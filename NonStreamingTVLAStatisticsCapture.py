# ============================================================
# STREAMING FULL TVLA + REQUIRED HEALTH + RAW SUM PERSIST
# ============================================================

import random
import os, time, json, pathlib, pickle, shutil
import numpy as np
import chipwhisperer as cw
from tqdm import tqdm
from math import comb

cw.scope_logger.setLevel(cw.logging.ERROR)

# ---------------- USER CONFIG ----------------
OUTDIR = "/home/adam/Desktop/Records/ThirdOrder"
N_TOTAL = 1_000_000
RUN_ID = 1

STREAM_BITS = 8
CHUNK_SAVE = 10_000
HEALTH_WINDOW = 1000

RETRIES_PER_TRACE = 5
ACK_TIMEOUT_MS = 1500
ADC_TIMEOUT_S = 8.0
TARGET_CLK_HZ = 24_000_000
USE_24MHZ_TARGET = True

FORCE_BAUD = None  # set to 921600 if you want to force, else None

SAMPLES_PER_OP = {
    0x30: 45000,
    0x31: 30000,
    0x32: 15000,
    0x33: 500,
    0x34: 85000,
    0x35: 75000,
    0x70: 45000

}

RUNS = [
    #{"cmd": 0x30, "name": "theta"},
    #{"cmd": 0x30, "name": "theta_null"},
    #{"cmd": 0x31, "name": "rho"},
    #{"cmd": 0x31, "name": "rho_null"},
    #{"cmd": 0x32, "name": "pi"},
    #{"cmd": 0x32, "name": "pi_null"},
    #{"cmd": 0x33, "name": "iota"},
    #{"cmd": 0x33, "name": "iota_null"},
    #{"cmd": 0x34, "name": "chi_null"},
    {"cmd": 0x35, "name": "chi_old"},
    {"cmd": 0x34, "name": "chi"},
    #{"cmd": 0x35, "name": "chi_from_scd_null"},
    #{"cmd": 0x70, "name": "nop"},   
]


FOLLOWUP_0x22_CMDS = {0x34, 0x35, 0x70}
FOLLOWUP_0x21_CMDS = {0x33}

EPS = 1e-20
# --------------------------------------------

pathlib.Path(OUTDIR).mkdir(parents=True, exist_ok=True)

# ================= helpers =================

def safe_name(s: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in s)

def op_dir(opname, cmd, samples):
    return pathlib.Path(OUTDIR) / f"{safe_name(opname)}_{cmd:02x}_{samples}_{N_TOTAL}_run{RUN_ID}_STREAMFULL"

def write_json(path: pathlib.Path, obj):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2))
    tmp.replace(path)

def append_jsonl(path: pathlib.Path, obj):
    with path.open("a") as f:
        f.write(json.dumps(obj) + "\n")

def save_pickle_atomic(path: pathlib.Path, obj):
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("wb") as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    tmp.replace(path)

def load_pickle(path: pathlib.Path):
    if path.exists():
        with path.open("rb") as f:
            return pickle.load(f)
    return None

def max_abs_and_argmax(arr: np.ndarray):
    a = np.abs(arr)
    idx = int(np.nanargmax(a))
    return float(a[idx]), idx

# ================= math =================

def centered_moment(m, c, k):
    out = 0.0
    for j in range(k + 1):
        out += comb(k, j) * ((-c) ** (k - j)) * m[j]
    return out

def centered_moment_2k(m, c, k):
    kk = 2 * k
    out = 0.0
    for j in range(kk + 1):
        out += comb(kk, j) * ((-c) ** (kk - j)) * m[j]
    return out

def welch_t(E0, E1, V0, V1, n0, n1):
    V0 = np.clip(V0, EPS, None)
    V1 = np.clip(V1, EPS, None)
    return (E0 - E1) / np.sqrt(V0 / n0 + V1 / n1 + 1e-12)

# ================= ChipWhisperer setup (ONCE) =================

scope  = cw.scope()
target = cw.target(scope, cw.targets.SimpleSerial2)

scope.default_setup()

scope.adc.stream_mode = False
#scope.adc.bits_per_sample = int(STREAM_BITS)

scope.clock.clkgen_src = "system"
scope.clock.adc_mul = 1
scope.clock.reset_dcms()

scope.trace.capture.trigger_source = "firmware trigger"

def resync_to_clk(new_clk_hz: int):
    old_clk  = float(scope.clock.clkgen_freq) or 7.37e6
    old_baud = int(getattr(target, "baud", 0) or 0)

    scope.clock.clkgen_freq = float(new_clk_hz)

    
    scope.io.nrst = "low";  time.sleep(0.05)
    scope.io.nrst = "high"; time.sleep(0.05)

    if FORCE_BAUD is not None:
        try:
            target.baud = int(FORCE_BAUD)
            return
        except Exception:
            pass

    ratio  = float(new_clk_hz) / float(old_clk)
    scaled = int(round(old_baud * ratio)) if old_baud else 921600
    for guess in (scaled, 748800, 768000, 921600):
        try:
            target.baud = int(guess)
            break
        except Exception:
            pass

def adc_errors_str():
    try:
        return str(scope.adc.errors)
    except Exception as ex:
        return f"<cannot read errors: {ex}>"

def clear_adc_errors():
    # Read-once often clears sticky status
    try:
        _ = scope.adc.errors
    except Exception:
        pass
    # Some builds expose a clearer
    for attr in ("clear_errors", "clear", "reset"):
        fn = getattr(scope.adc, attr, None)
        if callable(fn):
            try:
                fn()
                return
            except Exception:
                pass

def apply_stream(samples: int):
    scope.adc.stream_mode = True
    scope.adc.bits_per_sample = int(STREAM_BITS)
    scope.adc.presamples = 0
    scope.adc.offset = 0
    scope.adc.samples = int(samples)
    scope.adc.timeout = float(ADC_TIMEOUT_S)

if USE_24MHZ_TARGET:
    resync_to_clk(TARGET_CLK_HZ)

print(
    f"clkgen={float(scope.clock.clkgen_freq)/1e6:.2f}MHz "
    f"adc={float(scope.clock.adc_freq)/1e6:.2f}MHz "
    f"baud={int(getattr(target,'baud',-1))} "
    f"trigger={scope.trace.capture.trigger_source}"
)

# ================= firmware protocol =================

def prepare_state(pattern: int, cmd: int) -> float:
    t0 = time.perf_counter()

    target.flush()
    target.send_cmd(0x20, 0x00, bytes([pattern]))
    ack = target.read_cmd("e", timeout=ACK_TIMEOUT_MS)
    if not ack:
        raise RuntimeError("0x20 no ack")
    target.flush()

    if cmd in FOLLOWUP_0x22_CMDS:
        target.send_cmd(0x22, 0x00, b"\x00")
        ack2 = target.read_cmd("e", timeout=ACK_TIMEOUT_MS)
        if not ack2:
            raise RuntimeError("0x22 no ack")
        target.flush()

    if cmd in FOLLOWUP_0x21_CMDS:
        target.send_cmd(0x21, 0x00, b"\x00")
        ack2 = target.read_cmd("e", timeout=ACK_TIMEOUT_MS)
        if not ack2:
            raise RuntimeError("0x21 no ack")
        target.flush()

    return time.perf_counter() - t0

def capture_one(cmd: int, payload: bytes, samples: int):
    t0 = time.perf_counter()

    clear_adc_errors()
    target.flush()

    scope.arm()
    time.sleep(0.001)
    target.send_cmd(cmd, 0, payload)

    if scope.capture(poll_done=True):
        raise TimeoutError("no trigger")

    tr = scope.get_last_trace()
    if tr is None or len(tr) != samples:
        raise RuntimeError(f"trace_length_mismatch: got {0 if tr is None else len(tr)} expected {samples}")

    return np.asarray(tr, dtype=np.float64), (time.perf_counter() - t0)

# ================= main =================

for run in RUNS:
    cmd = run["cmd"]
    opname = run["name"]
    samples = SAMPLES_PER_OP[cmd]
    outdir = op_dir(opname, cmd, samples)
    outdir.mkdir(parents=True, exist_ok=True)

    resume_path = outdir / "resume.pkl"
    health_path = outdir / "health.jsonl"
    progress_path = outdir / "progress.json"
    meta_path = outdir / "meta.json"

    #apply_stream(samples)
    scope.adc.stream_mode = False
    scope.adc.samples = samples
    scope.adc.timeout = ADC_TIMEOUT_S
    scope.adc.presamples = 0
    scope.adc.offset = 0

    if not meta_path.exists():
        write_json(meta_path, dict(
            opname=opname,
            cmd=int(cmd),
            samples=int(samples),
            n_total=int(N_TOTAL),
            run_id=int(RUN_ID),
            stream_mode=False,
            stream_bits="default", #int(STREAM_BITS) for stream,
            adc_timeout_s=float(scope.adc.timeout),
            clkgen_freq=float(scope.clock.clkgen_freq),
            adc_freq=float(scope.clock.adc_freq),
            baud=int(getattr(target, "baud", -1)),
            trigger_source=str(scope.trace.capture.trigger_source),
            followup_0x22_cmds=sorted(list(FOLLOWUP_0x22_CMDS)),
            followup_0x21_cmds=sorted(list(FOLLOWUP_0x21_CMDS)),
            timestamp_start=time.time(),
        ))

    state = load_pickle(resume_path)
    if state is None:
        sumsF = [np.zeros(samples, dtype=np.float64) for _ in range(9)]
        sumsR = [np.zeros(samples, dtype=np.float64) for _ in range(9)]
        nF = nR = 0
        i = 0

        wall_start = time.time()
        capture_times = []
        ack_times = []

        global_mean_series = []
        global_var_series  = []

        energy_win = []
        mean_win = []
        var_win = []
        win_t0 = time.time()
        win_ok = 0

        e_sum_F = e_sum2_F = 0.0
        e_sum_R = e_sum2_R = 0.0
        e_min_F = float("inf"); e_max_F = float("-inf")
        e_min_R = float("inf"); e_max_R = float("-inf")

        counters = dict(
            retries_total=0,
            timeouts=0,
            prepare_failures=0,
            capture_failures=0,
            trace_length_mismatch=0,
            stream_errors=0,
        )
    else:
        (sumsF, sumsR, nF, nR, i,
         wall_start, capture_times, ack_times,
         global_mean_series, global_var_series,
         energy_win, mean_win, var_win,
         win_t0, win_ok,
         e_sum_F, e_sum2_F, e_sum_R, e_sum2_R,
         e_min_F, e_max_F, e_min_R, e_max_R,
         counters) = state

    print(f"\n=== {opname} cmd=0x{cmd:02x} samples={samples} (resume i={i}) ===")
    pbar = tqdm(total=N_TOTAL, initial=i)

    while i < N_TOTAL:
        group = random.getrandbits(1)
        payload = b"\x00"*16
        if "null" in opname:
            pattern = 0
            
        else:
            if group == 0:
                pattern = 1
            else:
                pattern = 0

        ok = False
        for attempt in range(RETRIES_PER_TRACE):
            try:
                ack_dt = prepare_state(pattern, cmd)

                tr, cap_dt = capture_one(cmd, payload, samples)

                # detect streaming errors (do not treat nonzero as hard-fail unless present)
                adc_err = getattr(scope.adc, "errors", None)
                try:
                    adc_nonzero = (int(adc_err) != 0)
                except Exception:
                    adc_nonzero = bool(adc_err)
                if adc_nonzero:
                    counters["stream_errors"] += 1

                ok = True
                break

            except TimeoutError:
                counters["timeouts"] += 1
                counters["capture_failures"] += 1
            except RuntimeError as e:
                s = str(e)
                if "trace_length_mismatch" in s:
                    counters["trace_length_mismatch"] += 1
                    counters["capture_failures"] += 1
                elif "no ack" in s:
                    counters["prepare_failures"] += 1
                else:
                    counters["capture_failures"] += 1
            except Exception:
                counters["capture_failures"] += 1
            finally:
                if not ok and attempt > 0:
                    counters["retries_total"] += 1
                if not ok and attempt < RETRIES_PER_TRACE - 1:
                    try:
                        target.flush()
                    except Exception:
                        pass
                    time.sleep(0.02)

        if not ok:
            save_pickle_atomic(resume_path, (
                sumsF, sumsR, nF, nR, i,
                wall_start, capture_times, ack_times,
                global_mean_series, global_var_series,
                energy_win, mean_win, var_win,
                win_t0, win_ok,
                e_sum_F, e_sum2_F, e_sum_R, e_sum2_R,
                e_min_F, e_max_F, e_min_R, e_max_R,
                counters
            ))
            write_json(progress_path, dict(
                opname=opname, cmd=int(cmd), samples=int(samples),
                next_index=int(i),
                n_fixed=int(nF), n_random=int(nR),
                aborted=True,
                adc_errors_last=adc_errors_str(),
                counters=counters,
                timestamp=time.time(),
            ))
            raise RuntimeError(f"Aborting: could not acquire trace at i={i} after {RETRIES_PER_TRACE} retries")

        # success
        ack_times.append(ack_dt)
        capture_times.append(cap_dt)

        tr_mean = float(tr.mean())
        tr_var  = float(tr.var())
        global_mean_series.append(tr_mean)
        global_var_series.append(tr_var)

        mean_win.append(tr_mean)
        var_win.append(tr_var)
        if len(mean_win) > HEALTH_WINDOW:
            mean_win.pop(0); var_win.pop(0)

        energy = float(np.dot(tr, tr))
        energy_win.append(energy)
        if len(energy_win) > HEALTH_WINDOW:
            energy_win.pop(0)

        if group == 0:
            e_sum_F += energy; e_sum2_F += energy * energy
            e_min_F = min(e_min_F, energy); e_max_F = max(e_max_F, energy)
        else:
            e_sum_R += energy; e_sum2_R += energy * energy
            e_min_R = min(e_min_R, energy); e_max_R = max(e_max_R, energy)

        # accumulate raw powers 1..8 
        p = np.ones(samples, dtype=np.float64)
        for k in range(1, 9):
            p *= tr
            if group == 0:
                sumsF[k] += p
            else:
                sumsR[k] += p

        if group == 0: nF += 1
        else:          nR += 1

        # window health log
        win_ok += 1
        if (i + 1) % HEALTH_WINDOW == 0:
            now = time.time()
            dtw = max(now - win_t0, 1e-9)
            recent_tps = win_ok / dtw

            append_jsonl(health_path, dict(
                i=int(i + 1),
                timestamp=now,
                n_fixed=int(nF),
                n_random=int(nR),
                recent_traces_per_second=float(recent_tps),
                recent_retry_rate=float(counters["retries_total"]) / float(max(nF + nR, 1)),
                window_energy_mean=float(np.mean(energy_win)) if energy_win else None,
                window_energy_std=float(np.std(energy_win)) if energy_win else None,
                window_global_mean=float(np.mean(mean_win)) if mean_win else None,
                window_global_var=float(np.mean(var_win)) if var_win else None,
                adc_errors_last=adc_errors_str(),
                counters=counters,
            ))
            win_t0 = now
            win_ok = 0

        # periodic persist
        if (i + 1) % CHUNK_SAVE == 0:
            save_pickle_atomic(resume_path, (
                sumsF, sumsR, nF, nR, i + 1,
                wall_start, capture_times, ack_times,
                global_mean_series, global_var_series,
                energy_win, mean_win, var_win,
                win_t0, win_ok,
                e_sum_F, e_sum2_F, e_sum_R, e_sum2_R,
                e_min_F, e_max_F, e_min_R, e_max_R,
                counters
            ))
            write_json(progress_path, dict(
                opname=opname, cmd=int(cmd), samples=int(samples),
                next_index=int(i + 1),
                n_fixed=int(nF), n_random=int(nR),
                adc_errors_last=adc_errors_str(),
                counters=counters,
                timestamp=time.time(),
            ))
            # optional backup like your working script
            if (i + 1) % 100_000 == 0:
                try:
                    shutil.copy(str(resume_path), str(resume_path) + ".bak")
                except Exception:
                    pass

        i += 1
        pbar.update(1)

    pbar.close()
    wall_end = time.time()

    # ================= Final stats =================
    if nF < 2 or nR < 2:
        raise RuntimeError(f"Need >=2 per group. Got nF={nF}, nR={nR}")

    mF = {0: 1.0, **{k: sumsF[k] / nF for k in range(1, 9)}}
    mR = {0: 1.0, **{k: sumsR[k] / nR for k in range(1, 9)}}

    mean_fixed  = mF[1]
    mean_random = mR[1]
    delta_mu = mean_fixed - mean_random

    var_fixed_raw  = (sumsF[2] - nF * mean_fixed**2) / (nF - 1)
    var_random_raw = (sumsR[2] - nR * mean_random**2) / (nR - 1)
    var_fixed  = np.clip(var_fixed_raw, EPS, None)
    var_random = np.clip(var_random_raw, EPS, None)

    N = nF + nR
    sum1_total = sumsF[1] + sumsR[1]
    sum2_total = sumsF[2] + sumsR[2]
    mu_total = sum1_total / N
    pooled_var = np.clip((sum2_total - N * mu_total**2) / (N - 1), EPS, None)

    rms_fixed  = np.sqrt(np.clip(mF[2], 0.0, None))
    rms_random = np.sqrt(np.clip(mR[2], 0.0, None))

    mu2F = np.clip(centered_moment(mF, mean_fixed, 2), EPS, None)
    mu3F = centered_moment(mF, mean_fixed, 3)
    mu4F = centered_moment(mF, mean_fixed, 4)
    mu2R = np.clip(centered_moment(mR, mean_random, 2), EPS, None)
    mu3R = centered_moment(mR, mean_random, 3)
    mu4R = centered_moment(mR, mean_random, 4)

    skew_fixed  = mu3F / (mu2F ** 1.5)
    skew_random = mu3R / (mu2R ** 1.5)
    kurt_fixed  = mu4F / (mu2F ** 2)
    kurt_random = mu4R / (mu2R ** 2)

    t1 = welch_t(mean_fixed, mean_random, var_fixed, var_random, nF, nR)

    t_var = (var_fixed - var_random) / np.sqrt(
        (2.0 * (var_fixed**2)) / (nF - 1) +
        (2.0 * (var_random**2)) / (nR - 1) + 1e-12
    )

    mu_pooled = (sumsF[1] + sumsR[1]) / N

    def save_arr(name, arr):
        np.save(outdir / f"{name}.npy", arr.astype(np.float64))

    save_arr("mean_fixed", mean_fixed)
    save_arr("mean_random", mean_random)
    save_arr("delta_mu", delta_mu)
    save_arr("var_fixed", var_fixed)
    save_arr("var_random", var_random)
    save_arr("var_fixed_raw", var_fixed_raw)
    save_arr("var_random_raw", var_random_raw)
    save_arr("pooled_var", pooled_var)
    save_arr("rms_fixed", rms_fixed)
    save_arr("rms_random", rms_random)
    save_arr("skew_fixed", skew_fixed)
    save_arr("skew_random", skew_random)
    save_arr("kurt_fixed", kurt_fixed)
    save_arr("kurt_random", kurt_random)
    save_arr("t1", t1)
    save_arr("t_var", t_var)

    t_peaks = {}
    m, ix = max_abs_and_argmax(t1)
    t_peaks["t1"] = dict(max_abs=m, index=ix)
    m, ix = max_abs_and_argmax(t_var)
    t_peaks["t_var"] = dict(max_abs=m, index=ix)

    for k in (2, 3, 4):
        # group-centered
        E0  = centered_moment(mF, mean_fixed, k)
        E1  = centered_moment(mR, mean_random, k)
        E20 = centered_moment_2k(mF, mean_fixed, k)
        E21 = centered_moment_2k(mR, mean_random, k)
        V0  = np.clip(E20 - E0**2, EPS, None)
        V1  = np.clip(E21 - E1**2, EPS, None)
        tkg = welch_t(E0, E1, V0, V1, nF, nR)
        save_arr(f"t{k}_group", tkg)
        m, ix = max_abs_and_argmax(tkg)
        t_peaks[f"t{k}_group"] = dict(max_abs=m, index=ix)

        # pooled-centered
        E0p  = centered_moment(mF, mu_pooled, k)
        E1p  = centered_moment(mR, mu_pooled, k)
        E20p = centered_moment_2k(mF, mu_pooled, k)
        E21p = centered_moment_2k(mR, mu_pooled, k)
        V0p  = np.clip(E20p - E0p**2, EPS, None)
        V1p  = np.clip(E21p - E1p**2, EPS, None)
        tkp = welch_t(E0p, E1p, V0p, V1p, nF, nR)
        save_arr(f"t{k}", tkp)
        m, ix = max_abs_and_argmax(tkp)
        t_peaks[f"t{k}"] = dict(max_abs=m, index=ix)

    write_json(outdir / "t_peaks.json", t_peaks)

    # energy stats
    def energy_stats(es, es2, n, emin, emax):
        if n <= 0:
            return dict(n=0)
        meanE = es / n
        varE  = max(es2 / n - meanE * meanE, 0.0)
        return dict(n=int(n), mean=float(meanE), std=float(np.sqrt(varE)), min=float(emin), max=float(emax))

    write_json(outdir / "trace_energy_stats.json", dict(
        energy_definition="sum(trace^2)",
        fixed=energy_stats(e_sum_F, e_sum2_F, nF, e_min_F, e_max_F),
        random=energy_stats(e_sum_R, e_sum2_R, nR, e_min_R, e_max_R),
    ))

    # timing stats
    cap = np.asarray(capture_times, dtype=np.float64)
    ack = np.asarray(ack_times, dtype=np.float64)

    np.save(outdir / "global_mean_series.npy", np.asarray(global_mean_series, dtype=np.float64))
    np.save(outdir / "global_var_series.npy",  np.asarray(global_var_series,  dtype=np.float64))

    write_json(outdir / "timing_stats.json", dict(
        wall_time_start=float(wall_start),
        wall_time_end=float(wall_end),
        total_wall_time=float(wall_end - wall_start),
        traces_per_second=float((nF + nR) / max(wall_end - wall_start, 1e-9)),
        mean_capture_time=float(np.mean(cap)) if cap.size else None,
        p95_capture_time=float(np.percentile(cap, 95)) if cap.size else None,
        mean_ack_time=float(np.mean(ack)) if ack.size else None,
        p95_ack_time=float(np.percentile(ack, 95)) if ack.size else None,
        counters=counters,
        n_fixed=int(nF),
        n_random=int(nR),
        adc_errors_last=adc_errors_str(),
    ))

    # persist raw sums
    save_pickle_atomic(outdir / "final_raw_sums.pkl", dict(
        sumsF=sumsF,
        sumsR=sumsR,
        nF=int(nF),
        nR=int(nR),
    ))

    # clean resume on success
    if resume_path.exists():
        resume_path.unlink()

print("\nComplete.")
