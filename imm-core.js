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
const EVW = 120 * FS;   // 흔들림 기록용 최근 120초 (50 Hz)
const EVB = [new Float32Array(EVW), new Float32Array(EVW), new Float32Array(EVW)];
let evi = 0;
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
    for (let c = 0; c < 3; c++) { BUF[c][wi] = o[c]; EVB[c][evi] = o[c]; } wi = (wi + 1) % W; evi = (evi + 1) % EVW; total++; lastOut = o; binSum = [0, 0, 0]; binN = 0; nextT += 1 / FS; }
  for (let c = 0; c < 3; c++) binSum[c] += lp2[c]; binN++;
}
/* I 안의 5단계: 상한 0.05, 0.1, 0.2, 0.4, 0.69 gal. 이름은 모두 I, 색은 무채색 명도 100·85·70·55·40%.
   소수 진도는 0.50~1.49를 0.2씩 나누고, 1단계(0~0.05 gal)는 직선, 2~5단계는 로그로 채운다. */
const I_SUB = [0.05, 0.1, 0.2, 0.4];   // 다섯 번째 상한은 II 하한(LOW[1])
const I_COLORS = ["#ffffff", "#d9d9d9", "#b3b3b3", "#8c8c8c", "#666666"];
function levelFrac(v) {
  if (!(v > 0)) return 0.5;
  let idx = 0; for (let i = 1; i < LOW.length; i++) if (v >= LOW[i]) idx = i;
  if (idx === 0) {
    if (v < I_SUB[0]) return 0.5 + 0.2 * v / I_SUB[0];
    const ub = [...I_SUB, LOW[1]];
    let j = 1; while (j < 4 && v >= ub[j]) j++;
    return 0.5 + 0.2 * j + 0.2 * Math.log(v / ub[j - 1]) / Math.log(ub[j] / ub[j - 1]);
  }
  const n = LOW.length - 1;
  if (idx >= n) return Math.min(n + 0.5 + Math.log(v / LOW[n]) / Math.log(LOW[n] / LOW[n - 1]), n + 1.49);
  return idx + 0.5 + Math.log(v / LOW[idx]) / Math.log(LOW[idx + 1] / LOW[idx]);
}
function classOf(f) { return Math.max(0, Math.min(9, Math.round(f) - 1)); }
function iSubOf(f) { return Math.max(1, Math.min(5, Math.floor((f - 0.5) / 0.2 + 1e-9) + 1)); }   // I 안의 단계 1~5
// 소수 진도 → 칩 색과 글자색
function colorOf(f) {
  const c = classOf(f);
  if (c === 0) { const s = iSubOf(f); return { bg: I_COLORS[s - 1], light: s === 5 }; }
  return { bg: MMI_COLORS[c], light: c >= 8 };
}
// 공유 코드: 0 = 잡음 이하, 1~5 = I의 1~5단계, 6~14 = II~X+
function codeOf(f) { if (!(f > 0)) return 0; const c = classOf(f); return c === 0 ? iSubOf(f) : c + 5; }
function codeInfo(k) {
  if (!(k > 0)) return null;
  if (k <= 5) return { name: "I", bg: I_COLORS[k - 1], light: k === 5, sub: k };
  const c = Math.min(9, k - 5); return { name: CLASSES[c], bg: MMI_COLORS[c], light: c >= 8 };
}
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

/* 건물 보정: 균일 전단 기둥 1~5차 모드로 지금 층의 수평 증폭 H(f)를 계산해 나눈다 (역보정 최대 3배) */
const bld = { on: false, top: 25, cur: 20, T1: null, zeta: 5 };
let HinvR = null, HinvI = null;
function estT1(top) { return 0.0488 * Math.pow(2.8 * top, 0.75); }   // 국내 벽식 아파트 기준식 (층고 2.8 m)
function bldT1() { return bld.T1 || estT1(bld.top); }
function bldH(fa) {
  const T1 = bldT1(), z = bld.zeta / 100, zr = Math.max(0, Math.min(1, (bld.cur - 1) / bld.top));
  let hr = 1, hi = 0;
  if (fa > 0) for (let n = 1; n <= 5; n++) {
    const m = 2 * n - 1, r = fa * T1 / m, G = 4 * (n % 2 ? 1 : -1) / (m * Math.PI), ph = Math.sin(m * Math.PI * zr / 2);
    const a = 1 - r * r, b = 2 * z * r, d2 = a * a + b * b, c = G * ph * r * r / d2;
    hr += c * a; hi -= c * b;
  }
  return [hr, hi];
}
function buildHinv() {
  if (!bld.on) { HinvR = HinvI = null; return; }
  HinvR = new Float64Array(N); HinvI = new Float64Array(N);
  for (let k = 0; k < N; k++) {
    const sf = (k <= N / 2 ? k : k - N) * FS / N;
    let [hr, hi] = bldH(Math.abs(sf)); if (sf < 0) hi = -hi;
    const d = hr * hr + hi * hi; let ir = hr / d, ii = -hi / d;
    const mg = Math.sqrt(ir * ir + ii * ii); if (mg > 3) { ir *= 3 / mg; ii *= 3 / mg; }
    HinvR[k] = ir; HinvI[k] = ii;
  }
}
// 조용할 때 수평 미세진동 스펙트럼의 0.15~3 Hz 봉우리 → 건물 고유주기 (20초 이상 데이터 필요)
function measureT1() {
  const Lw = Math.min(total, W); if (Lw < 20 * FS) return null;
  const re = new Float64Array(N), im = new Float64Array(N);
  for (let c = 0; c < 2; c++) {
    const a = new Float64Array(Lw); for (let k = 0; k < Lw; k++) a[k] = BUF[c][(wi - Lw + k + W) % W];
    let st = 0, sx = 0, stt = 0, stx = 0; for (let k = 0; k < Lw; k++) { st += k; sx += a[k]; stt += k * k; stx += k * a[k]; }
    const sl = (Lw * stx - st * sx) / (Lw * stt - st * st), ic = (sx - sl * st) / Lw;
    for (let k = 0; k < Lw; k++) (c ? im : re)[k] = a[k] - ic - sl * k;
  }
  fft(re, im, false);
  const P = new Float64Array(N / 2);
  for (let k = 1; k < N / 2; k++) P[k] = (re[k] ** 2 + im[k] ** 2 + re[N - k] ** 2 + im[N - k] ** 2) / 2;
  const S = new Float64Array(N / 2);
  for (let k = 3; k < N / 2 - 3; k++) { let q = 0; for (let j = -2; j <= 2; j++) q += P[k + j]; S[k] = q / 5; }
  const k0 = Math.ceil(0.15 * N / FS), k1 = Math.floor(3 * N / FS);
  let kb = k0; for (let k = k0; k <= k1; k++) if (S[k] > S[kb]) kb = k;
  const band = Array.from(S.slice(k0, k1 + 1)).sort((a, b) => a - b), med = band[band.length >> 1];
  if (!(S[kb] > 4 * med) || kb === k0 || kb === k1) return null;
  const y0 = S[kb - 1], y1 = S[kb], y2 = S[kb + 1], dk = (y0 - y2) / (2 * (y0 - 2 * y1 + y2)) || 0;
  return { T: 1 / ((kb + dk) * FS / N), ratio: S[kb] / med };
}
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
  if (HinvR) for (let k = 0; k < N; k++) {   // 수평(x+iy)만 건물 증폭을 나눈다
    const r = R1[k] * HinvR[k] - I1[k] * HinvI[k], i = R1[k] * HinvI[k] + I1[k] * HinvR[k]; R1[k] = r; I1[k] = i;
  }
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

/* ===== 흔들림 기록 파형 압축 =====
   최근 n개 표본(50 Hz)을 축마다 b칸으로 나눠 칸별 최솟값·최댓값만 남기고, μ-law 8비트로 줄인 뒤 deflate로 압축한다.
   형식 v1: { v: 1, b: 칸 수, dur: 초, s: [축별 최대 |가속도| gal], z: 1(deflate) 또는 0, d: base64 }
   바이트 순서: 축 X, Y, Z 각각 칸마다 [최솟값, 최댓값] (0~254, 127이 0) */
const MU = 255;
function muEnc(x) { const y = Math.sign(x) * Math.log1p(MU * Math.min(1, Math.abs(x))) / Math.log1p(MU); return Math.round(y * 127) + 127; }
function muDec(q) { const y = (q - 127) / 127; return Math.sign(y) * (Math.pow(1 + MU, Math.abs(y)) - 1) / MU; }
function recentRaw(n) {
  n = Math.min(n, total, EVW);
  return [0, 1, 2].map(c => { const a = new Float64Array(n); for (let k = 0; k < n; k++) a[k] = EVB[c][(evi - n + k + EVW) % EVW]; return a; });
}
function b64(bytes) { let s = ""; for (let i = 0; i < bytes.length; i += 0x8000) s += String.fromCharCode.apply(null, bytes.subarray(i, i + 0x8000)); return btoa(s); }
function unb64(str) { const s = atob(str), a = new Uint8Array(s.length); for (let i = 0; i < s.length; i++) a[i] = s.charCodeAt(i); return a; }
async function pipeBytes(bytes, stream) { const r = new Blob([bytes]).stream().pipeThrough(stream); return new Uint8Array(await new Response(r).arrayBuffer()); }
async function encodeWave(n, bins) {
  const xs = recentRaw(n), len = xs[0].length, b = Math.min(bins || 600, len);
  const out = new Uint8Array(6 * b), sc = [];
  xs.forEach((a, c) => {
    let st = 0, sx = 0, stt = 0, stx = 0; for (let k = 0; k < len; k++) { st += k; sx += a[k]; stt += k * k; stx += k * a[k]; }
    const sl = (len * stx - st * sx) / (len * stt - st * st), ic = (sx - sl * st) / len;
    let pk = 0; for (let k = 0; k < len; k++) { a[k] -= ic + sl * k; pk = Math.max(pk, Math.abs(a[k])); }   // 중력·기울기 성분 제거
    pk = pk || 1e-6; sc.push(Number(pk.toPrecision(4)));
    for (let j = 0; j < b; j++) {
      const k0 = Math.floor(j * len / b), k1 = Math.max(k0 + 1, Math.floor((j + 1) * len / b));
      let mn = Infinity, mx = -Infinity; for (let k = k0; k < k1; k++) { if (a[k] < mn) mn = a[k]; if (a[k] > mx) mx = a[k]; }
      out[(c * b + j) * 2] = muEnc(mn / pk); out[(c * b + j) * 2 + 1] = muEnc(mx / pk);
    }
  });
  let z = 0, data = out;
  if (typeof CompressionStream !== "undefined") { try { data = await pipeBytes(out, new CompressionStream("deflate")); z = 1; } catch (e) { data = out; } }
  return { v: 1, b, dur: Math.round(len / FS * 10) / 10, s: sc, z, d: b64(data) };
}
async function decodeWave(w) {
  let bytes = unb64(w.d);
  if (w.z) bytes = await pipeBytes(bytes, new DecompressionStream("deflate"));
  return [0, 1, 2].map(c => {
    const mn = new Float32Array(w.b), mx = new Float32Array(w.b);
    for (let j = 0; j < w.b; j++) { mn[j] = muDec(bytes[(c * w.b + j) * 2]) * w.s[c]; mx[j] = muDec(bytes[(c * w.b + j) * 2 + 1]) * w.s[c]; }
    return { mn, mx, peak: w.s[c] };
  });
}
