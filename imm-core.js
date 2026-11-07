/* IMM 계산 핵심 (IMM 진도계와 같은 방식)
   pushRaw(t초, x, y, z gal)로 넣고 1초마다 computeIMM(잡음바닥)을 부른다.
   반환 v = [IMM_(0.2,0.5), IMM_(0.5,1.5), IMM_(1.5,7.5), 대표 IMM] 소수 진도 (0 = 잡음 이하), T = 대표 주기 */
/* ===================== IMM 계산 (IMM 진도계와 같은 방식) ===================== */
const FS = 50, W = 2048, N = 4096, T0 = 1.616, K = 4, TMIN = 0.1, TMAX = 10, PER_OCT = 6;
const CLASSES = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X+"];
const LOW = [0, 0.07, 0.23, 0.76, 2.56, 6.86, 14.73, 31.66, 68.01, 146.14].map(p => p / 100 * 980.665);
const MMI_COLORS = ["#ffffff", "#bfccff", "#a0e6ff", "#80ffff", "#7aff93", "#ffff00", "#ffc800", "#ff9100", "#ff0000", "#c80000"];
const NT = Math.round(Math.log2(TMAX / TMIN) * PER_OCT) + 1;
const TS = Array.from({ length: NT }, (_, i) => TMIN * Math.pow(2, i / PER_OCT));
const RANGES = [{ a: 0.2, b: 0.5, k: "IMM_(0.2,0.5)" }, { a: 0.5, b: 1.5, k: "IMM_(0.5,1.5)" }, { a: 1.5, b: 7.5, k: "IMM_(1.5,7.5)" }];
const ENDS = [0.2, 0.5, 1.5, 7.5];

function hcut(f) { const x2 = (f / 10) * (f / 10); return Math.pow(1 + x2 * (0.694 + x2 * (0.241 + x2 * (0.0557 + x2 * (0.009664 + x2 * (0.00134 + x2 * 0.000155))))), -0.5); }
function jmaGainF(f) { return f <= 0 ? 0 : Math.sqrt(1 / f) * hcut(f) * Math.sqrt(1 - Math.exp(-Math.pow(f / 0.5, 3))); }
const FMAX = jmaGainF(1 / T0);
function nfac(T) { const w = 1 / (1 + Math.pow(T0 / T, K)); return Math.pow(T0 / T, w); }
function sigWeight(f) { return 1 / nfac(Math.min(1 / Math.max(f, 1e-12), TMAX)); }
function binFreq(k) { return (k <= N / 2 ? k : N - k) * FS / N; }
function immGain(T) { const g = new Float64Array(N); for (let k = 0; k < N; k++) { const f = binFreq(k); g[k] = f > 0 ? jmaGainF(f * T / T0) / FMAX * sigWeight(f) : 0; } return g; }
const GAINS = TS.map(immGain), EGAINS = ENDS.map(immGain);

const BITREV = new Uint32Array(N), COS = new Float64Array(N / 2), SIN = new Float64Array(N / 2);
(function () { const bits = Math.log2(N); for (let i = 0; i < N; i++) { let r = 0; for (let b = 0; b < bits; b++) r = (r << 1) | ((i >> b) & 1); BITREV[i] = r; }
  for (let i = 0; i < N / 2; i++) { COS[i] = Math.cos(2 * Math.PI * i / N); SIN[i] = Math.sin(2 * Math.PI * i / N); } })();
function fft(re, im, inv) {
  for (let i = 0; i < N; i++) { const j = BITREV[i]; if (j > i) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; } }
  const sg = inv ? 1 : -1;
  for (let size = 2; size <= N; size <<= 1) { const h = size >> 1, st = N / size;
    for (let s0 = 0; s0 < N; s0 += size) for (let j = 0; j < h; j++) { const c = COS[j * st], s = sg * SIN[j * st], a = s0 + j, b = a + h;
      const tr = re[b] * c - im[b] * s, ti = re[b] * s + im[b] * c; re[b] = re[a] - tr; im[b] = im[a] - ti; re[a] += tr; im[a] += ti; } }
  if (inv) for (let i = 0; i < N; i++) { re[i] /= N; im[i] /= N; }
}

const BUF = [new Float32Array(W), new Float32Array(W), new Float32Array(W)];
let wi = 0, total = 0, nextT = null, lastRawT = 0, binSum = [0, 0, 0], binN = 0, lastOut = [0, 0, 0], lp1 = [0, 0, 0], lp2 = [0, 0, 0];
const LP_TAU = 1 / (2 * Math.PI * 25);
function resetBuf() { wi = 0; total = 0; nextT = null; lastRawT = 0; binN = 0; }
function pushRaw(t, x, y, z) {   // 25 Hz 저역통과 후 20 ms 구간 평균으로 50 Hz 재샘플링
  const v = [x, y, z];
  if (nextT === null || t - lastRawT > 1 || t < lastRawT) { nextT = t + 1 / FS; binSum = v.slice(); binN = 1; lastRawT = t; lastOut = v.slice(); lp1 = v.slice(); lp2 = v.slice(); return; }
  const al = (t - lastRawT) / (LP_TAU + (t - lastRawT));
  for (let c = 0; c < 3; c++) { lp1[c] += al * (v[c] - lp1[c]); lp2[c] += al * (lp1[c] - lp2[c]); }
  lastRawT = t;
  while (t >= nextT) { const o = binN ? binSum.map(q => q / binN) : lastOut;
    for (let c = 0; c < 3; c++) BUF[c][wi] = o[c]; wi = (wi + 1) % W; total++; lastOut = o; binSum = [0, 0, 0]; binN = 0; nextT += 1 / FS; }
  for (let c = 0; c < 3; c++) binSum[c] += lp2[c]; binN++;
}
function levelFrac(v) {
  if (!(v > 0)) return 0.5;
  let idx = 0; for (let i = 1; i < LOW.length; i++) if (v >= LOW[i]) idx = i;
  if (idx === 0) return 0.5 + Math.min(v / LOW[1], 1) * 0.999;
  const n = LOW.length - 1;
  if (idx >= n) return Math.min(n + 0.5 + Math.log(v / LOW[n]) / Math.log(LOW[n] / LOW[n - 1]), n + 1.49);
  return idx + 0.5 + Math.log(v / LOW[idx]) / Math.log(LOW[idx + 1] / LOW[idx]);
}
function classOf(f) { return Math.max(0, Math.min(9, Math.round(f) - 1)); }
function fmtFrac(v) { return (Math.floor(v * 100 + 1e-9) / 100).toFixed(2); }

const R1 = new Float64Array(N), I1 = new Float64Array(N), R2 = new Float64Array(N), I2 = new Float64Array(N);
const Ar = new Float64Array(N), Ai = new Float64Array(N), Br = new Float64Array(N), Bi = new Float64Array(N);
function peakFor(g, Lw, from) {
  for (let k = 0; k < N; k++) { Ar[k] = R1[k] * g[k]; Ai[k] = I1[k] * g[k]; Br[k] = R2[k] * g[k]; Bi[k] = I2[k] * g[k]; }
  fft(Ar, Ai, true); fft(Br, Bi, true);
  let p = 0; for (let i = from; i < Lw; i++) { const m = Ar[i] * Ar[i] + Ai[i] * Ai[i] + Br[i] * Br[i]; if (m > p) p = m; }
  return Math.sqrt(p);
}
// 최근 40초 창 → [IMM_(0.2,0.5), IMM_(0.5,1.5), IMM_(1.5,7.5), 대표 IMM] 소수 진도와 대표 주기 T*
// noise: 주기별 잡음 바닥(TS와 같은 길이). 그 1.5배를 넘지 못한 주기는 판정에서 뺀다. 모두 빠지면 0 (잡음 이하)
const NOISE_MARGIN = 1.5;
function noiseAt(arr, T) {
  if (!arr) return null;
  if (T <= TS[0]) return arr[0];
  if (T >= TS[NT - 1]) return arr[NT - 1];
  const x = Math.log2(T / TMIN) * PER_OCT, i = Math.floor(x), r = x - i;
  return arr[i] * Math.pow(arr[i + 1] / arr[i], r);
}
function computeIMM(noise) {
  const Lw = Math.min(total, W);
  if (Lw < 5 * FS) return null;
  R1.fill(0); I1.fill(0); R2.fill(0); I2.fill(0);
  const xs = [0, 1, 2].map(c => { const a = new Float64Array(Lw); for (let k = 0; k < Lw; k++) a[k] = BUF[c][(wi - Lw + k + W) % W]; return a; });
  for (const a of xs) {
    let st = 0, sx = 0, stt = 0, stx = 0; for (let k = 0; k < Lw; k++) { st += k; sx += a[k]; stt += k * k; stx += k * a[k]; }
    const sl = (Lw * stx - st * sx) / (Lw * stt - st * st), ic = (sx - sl * st) / Lw; for (let k = 0; k < Lw; k++) a[k] -= ic + sl * k;
    const nt = Math.floor(Lw * 0.05); for (let k = 0; k < nt; k++) a[k] *= 0.5 - 0.5 * Math.cos(Math.PI * k / nt);
    const ne = Math.min(150, Math.floor(Lw / 4)); for (let k = 0; k < ne; k++) a[Lw - 1 - k] *= 0.5 - 0.5 * Math.cos(Math.PI * k / ne);
  }
  for (let k = 0; k < Lw; k++) { R1[k] = xs[0][k]; I1[k] = xs[1][k]; R2[k] = xs[2][k]; }
  fft(R1, I1, false); fft(R2, I2, false);
  const from = Math.floor(Lw * 0.05);
  const A = GAINS.map(g => peakFor(g, Lw, from)), E = EGAINS.map(g => peakFor(g, Lw, from));
  const ok = (a, n) => !(n !== null && n !== undefined && a < NOISE_MARGIN * n);
  const lv = A.map((a, i) => ok(a, noise ? noise[i] : null) ? levelFrac(a) : 0);
  const ev = E.map((a, j) => ok(a, noiseAt(noise, ENDS[j])) ? levelFrac(a) : 0);
  let best = 0; lv.forEach((v, i) => { if (v > lv[best]) best = i; });
  const ranges = RANGES.map(r => {
    let f = Math.max(ev[ENDS.indexOf(r.a)], ev[ENDS.indexOf(r.b)]);
    TS.forEach((T, i) => { if (T > r.a && T < r.b && lv[i] > f) f = lv[i]; });
    return f;
  });
  return { A, v: [...ranges, lv[best]], T: lv[best] > 0 ? TS[best] : 0 };
}
