# I-MAIHDA HIC-MIC Simulation v3.1

> **Tiếng Việt** | [English version below](#english)

Bộ công cụ mô phỏng dữ liệu tổng hợp kiểm định độ nhạy của **VPC** và **PCV** — hai chỉ số cốt lõi trong I-MAIHDA — trước tác động của tỉ lệ hiện mắc, strata thưa, và sai số phát hiện có khuôn mẫu SES. Repository gồm **Python** (bản chính) và **R package `imaihda`** (bản tái lập hoàn chỉnh, có thể cài đặt).

⚠️ **Không sử dụng dữ liệu thật.** Không có tuyên bố thực nghiệm về bất kỳ quần thể nào. Đây là minh chứng phương pháp luận.

---

##  Mục lục

1. [Câu hỏi nghiên cứu](#câu-hỏi-nghiên-cứu)
2. [Phương pháp](#phương-pháp)
3. [Kịch bản kiểm định](#kịch-bản-kiểm-định)
4. [Kết quả benchmark](#kết-quả-benchmark)
5. [Hình minh họa](#hình-minh-họa)
6. [R package `imaihda`](#r-package-imaihda)
7. [Đối chứng Python ↔ R](#đối-chứng-python--r)
8. [So sánh R package với R scripts cũ](#so-sánh-r-package-với-r-scripts-cũ)
9. [FAQs](#faqs)
10. [Tài liệu tham khảo](#tài-liệu-tham-khảo)

---

## Câu hỏi nghiên cứu

Nếu một cohort thu nhập trung bình (MIC) cho thấy **VPC cao hơn** hoặc **PCV thấp hơn** so với cohort thu nhập cao (HIC), liệu điều đó có nhất thiết nghĩa là cấu trúc bất bình đẳng giao thoa khác biệt?

**Trả lời: Không.** VPC và PCV có thể dao động theo tỉ lệ hiện mắc, strata thưa, và sai số phát hiện có khuôn mẫu SES — ngay cả khi cấu trúc giao thoa thực sự không đổi.

---

## Phương pháp

Mô phỏng cá thể lồng trong **36 strata giao thoa** (giới tính × học vấn × tài sản × nông thôn). Tính toán chẩn đoán nhanh qua **empirical logit** và **weighted method-of-moments**, không dùng GLMM đầy đủ.

### Công thức

**VPC** (Variance Partition Coefficient) trên thang latent logistic:

$$\text{VPC} = \frac{\sigma^2_{\text{stratum}}}{\sigma^2_{\text{stratum}} + \pi^2/3} \times 100\%$$

**PCV** (Proportional Change in Variance):

$$\text{PCV} = \frac{\sigma^2_{\text{null}} - \sigma^2_{\text{main}}}{\sigma^2_{\text{null}}} \times 100\%$$

Trong đó:
- `σ²_null` = phương sai giữa các strata từ mô hình null (chỉ có strata, không có biến giải thích)
- `σ²_main` = phương sai giữa các strata còn lại sau khi thêm hiệu ứng chính cộng-gộp
- `π²/3 ≈ 3.29` = phương sai mức cá thể của phân phối logistic chuẩn

---

## Kịch bản kiểm định

| Kịch bản | Mục đích |
|----------|----------|
| **A** | Gradient xã hội thuần cộng-gộp, phát hiện đồng đều |
| **B** | Có tương tác giao thoa thực sự, phát hiện đồng đều |
| **C** | Cấu trúc cộng-gộp với sai số phát hiện theo khuôn mẫu SES |
| **D** | Tương tác giao thoa + sai số phát hiện SES |
| **E** | Tương tác giao thoa, bệnh hiếm, strata thưa |

---

## Kết quả benchmark

### Python (PCG64 RNG)

| Scenario | Tỉ lệ hiện mắc | VPC null | VPC main | PCV | Strata nhỏ nhất |
|----------|:--------------:|:--------:|:--------:|:---:|:---------------:|
| A | 23.3% | 4.32 | 0.00 | **100.0** | 144 |
| B | 27.1% | 22.58 | 15.78 | **35.8** | 144 |
| C | 11.3% | 0.00 | 0.00 | NaN | 144 |
| D | 13.7% | 13.68 | 8.80 | **39.1** | 144 |
| E | 9.1% | 14.70 | 9.44 | **39.5** | 1 |

### R (`imaihda` package, Mersenne Twister RNG)

| Scenario | Tỉ lệ hiện mắc | VPC null | VPC main | PCV | Strata nhỏ nhất |
|----------|:--------------:|:--------:|:--------:|:---:|:---------------:|
| A | 23.6% | 4.50 | 0.00 | **100.0** | 130 |
| B | 26.3% | 17.20 | 9.12 | **51.7** | 130 |
| C | 11.6% | 0.82 | 0.18 | **78.0** | 130 |
| D | 13.0% | 16.06 | 8.02 | **54.5** | 130 |
| E | 11.6% | 19.47 | 2.98 | **87.3** | 2 |

### Pass/fail benchmarks — Cả Python và R đều PASS cả 6/6

| # | Tiêu chí | Python | R |
|---|----------|:------:|:--:|
| 1 | A thuần cộng-gộp: PCV ≥ 80, VPC_main < 1 | ✅ | ✅ |
| 2 | B có tương tác thực: VPC_null(B) > VPC_null(A) + 5 | ✅ | ✅ |
| 3 | B để lại phương sai dư: PCV < 70 | ✅ | ✅ |
| 4 | C sai số phát hiện làm giảm tỉ lệ hiện mắc | ✅ | ✅ |
| 5 | D sai số phát hiện che lấp VPC tương tác | ✅ | ✅ |
| 6 | E strata thưa được gắn cờ | ✅ | ✅ |

> **Kết luận chính:** Cả Python và R đều xác nhận rằng VPC/PCV nhạy với prevalence, strata thưa, và sai số phát hiện. Không thể so sánh HIC‑MIC một cách thô nếu không có chẩn đoán đi kèm.

---

## Hình minh họa

### 1. Phân bố VPC‑PCV theo kịch bản

![VPC-PCV scenario plot](figures/scenario_vpc_pcv.png)

**Diễn giải:** Kịch bản A (góc trên bên trái) thể hiện cấu trúc thuần cộng-gộp (PCV=100%). B dịch chuyển sang phải (VPC cao hơn) và xuống dưới (PCV thấp hơn) do có tương tác giao thoa thực sự. D cho thấy sai số phát hiện có thể che lấp VPC của tương tác. E minh họa tác động của strata thưa.

### 2. Quét sai số phát hiện (detection sweep)

![Detection sweep](figures/detection_sweep.png)

**Diễn giải:** Khi cường độ sai số phát hiện theo khuôn mẫu SES tăng lên, tỉ lệ hiện mắc quan sát giảm (đường đứt nét). VPC ban đầu giảm (che lấu) rồi sau đó có thể tăng trở lại ở mức sai số rất cao — cho thấy tác động không đơn điệu.

---

## R package `imaihda`

###  Cài đặt

```r
# Cài từ GitHub
remotes::install_github("nguyenminh2301/-i-maihda", subdir = "imaihda")

# Hoặc clone về cài local
# git clone https://github.com/nguyenminh2301/-i-maihda.git
# devtools::install("duong-dan/-i-maihda/imaihda")
```

**Yêu cầu:** R ≥ 4.0, packages: `stats` (base), `ggplot2` (gợi ý), `testthat` (gợi ý).

###  Load package

```r
library(imaihda)
```

###  Sử dụng từng hàm

#### `vpc_latent()` — Tính VPC từ phương sai strata

```r
vpc_latent(0.5)     # 13.19% — VPC khi σ²_stratum = 0.5
vpc_latent(0)       # 0%
vpc_latent(pi^2/3)  # 50% — VPC = 50% khi σ² = π²/3
```

#### `pcv()` — Tính PCV từ hai phương sai

```r
pcv(1.0, 0.25)  # 75% — hầu hết phương sai được giải thích bởi hiệu ứng chính
pcv(0.5, 0.4)   # 20% — còn nhiều tương tác dư
pcv(0, 0)       # NaN — không xác định khi phương sai null ≤ 0
```

#### `simulate_intersectional_data()` — Sinh dữ liệu tổng hợp

```r
# Kịch bản cơ bản (thuần cộng-gộp)
df <- simulate_intersectional_data(n = 2000, seed = 42)

# Có tương tác giao thoa
df_b <- simulate_intersectional_data(n = 2000, interaction_sd = 0.9, seed = 42)

# Có sai số phát hiện SES
df_c <- simulate_intersectional_data(n = 2000, detection_strength = 0.8, seed = 42)

# Strata thưa, bệnh hiếm
df_e <- simulate_intersectional_data(
  n = 1000, prevalence_shift = -3.0, interaction_sd = 0.9,
  sparse = TRUE, seed = 42
)

# So sánh tỉ lệ hiện mắc quan sát và thực
mean(df_c$y)       # quan sát (thấp hơn do sai số phát hiện)
mean(df_c$y_true)  # thực
```

#### `fit_imaihda()` — Chẩn đoán MAIHDA một lần gọi

```r
df <- simulate_intersectional_data(n = 3000, seed = 123)
res <- fit_imaihda(df)

res$n_strata              # 36
res$overall_prevalence    # ~0.23
res$vpc_null              # VPC từ mô hình null
res$vpc_main              # VPC sau hiệu ứng chính
res$pcv                   # PCV (%)
res$var_null              # Phương sai giữa strata (null)
res$var_main              # Phương sai giữa strata (main)
res$min_stratum_n         # Cỡ strata nhỏ nhất
```

#### `scenario_grid()` + `evaluate_benchmarks()` — Chạy toàn bộ 5 kịch bản

```r
grid <- scenario_grid()
results <- do.call(rbind, lapply(names(grid), function(nm) {
  res <- fit_scenario(nm, grid[[nm]])
  as.data.frame(res)
}))
benchmarks <- evaluate_benchmarks(results)
print(benchmarks)  # 6 dòng pass/fail
```

###  Chạy tests

```r
devtools::test("imaihda")   # 39 testthat tests
```

---

## Đối chứng Python ↔ R

### Khác biệt về RNG

| | Python | R |
|---|---|---|
| **Engine** | PCG64 (`numpy.random.default_rng`) | Mersenne Twister (`set.seed`) |
| **Hạt giống** | 42 | 42 |
| **Kết quả số** | Khác | Khác |
| **Kết luận benchmark** | Giống (6/6 pass) | Giống (6/6 pass) |

###  So sánh từng chỉ số

| Chỉ số | Python (điển hình) | R (điển hình) | Khớp định tính? |
|--------|:------------------:|:-------------:|:---------------:|
| VPC_null(A) | 4.32 | 4.50 | ✅ Thấp |
| VPC_null(B) | 22.58 | 17.20 | ✅ Cao hơn A |
| PCV(A) | 100.0 | 100.0 | ✅ Cộng-gộp hoàn toàn |
| PCV(B) | 35.8 | 51.7 | ✅ < 70 (còn tương tác dư) |
| Tỉ lệ C/A | 11.3/23.3 | 11.6/23.6 | ✅ Giảm ~50% |
| VPC(D) < VPC(B) | Có (13.68 < 22.58) | Có (16.06 < 17.20) | ✅ Che lấp |
| Min_n(E) | 1 | 2 | ✅ Thưa |

> Cả hai ngôn ngữ đều đi đến **cùng kết luận định tính**. Khác biệt về số liệu tuyệt đối là do RNG, không ảnh hưởng đến diễn giải khoa học.

---

## So sánh R package với R scripts cũ

| Tiêu chí | R scripts cũ (v3.1) | R package `imaihda` (v3.2) |
|----------|---------------------|---------------------------|
| **Cấu trúc** | File `.R` rời rạc trong `R/` | Package chuẩn với DESCRIPTION, NAMESPACE |
| **Cài đặt** | `source("R/simulate.R")` thủ công | `install_github()` hoặc `devtools::install()` |
| **Tài liệu** | Comment nội bộ | Roxygen2 với `@examples`, `@references` |
| **Hàm exported** | Không phân biệt | 6 hàm public, 6 hàm internal |
| **Tests** | 4 test_that blocks (thủ công) | 39 testthat assertions (tự động) |
| **Tương thích** | Gắn với thư mục dự án WZB | Độc lập, dùng được ở mọi dự án |
| **Kết quả** | Khớp với Python về benchmark | **Giống hệt** scripts cũ (cùng code, cùng RNG) |

> **Tính nhất quán:** R package dùng **cùng thuật toán** với R scripts cũ. Với cùng seed, kết quả **giống hệt 100%** vì logic tính toán không đổi — chỉ khác về cách tổ chức code.

---

## FAQs

<details>
<summary><b>1. Đây có phải là công cụ ước lượng mới cho MAIHDA không?</b></summary>

Không. Đây là **minh chứng phương pháp luận** dùng empirical-logit diagnostic nhanh để stress-test VPC/PCV. Trong nghiên cứu thực nghiệm, cần dùng GLMM đầy đủ (ví dụ: `lme4::glmer` trong R hoặc `statsmodels.MixedLM` trong Python).
</details>

<details>
<summary><b>2. Tại sao kết quả Python và R khác nhau?</b></summary>

Vì **RNG engine khác nhau** (PCG64 vs Mersenne Twister). Cùng seed `42` nhưng chuỗi số ngẫu nhiên khác → dữ liệu mô phỏng khác → số liệu VPC/PCV khác. Tuy nhiên, **kết luận định tính và benchmark pass/fail giống hệt nhau**. Đây là đặc điểm mong đợi của bất kỳ cross-language reproduction nào dùng RNG.
</details>

<details>
<summary><b>3. Tôi có thể dùng package này với dữ liệu thật không?</b></summary>

Có thể dùng `fit_imaihda()` để chẩn đoán nhanh, nhưng phải hiểu rằng đây là **xấp xỉ chẩn đoán** (method-of-moments), không phải ước lượng GLMM đầy đủ. Đối với dữ liệu thật và công bố khoa học, dùng `lme4::glmer(y ~ (1|stratum) + covariates, family=binomial)`.
</details>

<details>
<summary><b>4. Làm sao để sinh lại hình (figures)?</b></summary>

```r
library(imaihda)
library(ggplot2)

# Chạy tất cả kịch bản
grid <- scenario_grid()
rows <- lapply(names(grid), function(nm) fit_scenario(nm, grid[[nm]]))
results <- do.call(rbind, lapply(rows, as.data.frame))

# Vẽ VPC-PCV scatter
ggplot(results, aes(vpc_null, pcv, label = scenario)) +
  geom_point(size = 3, color = "#21918c") +
  geom_text(hjust = -0.3, family = "serif") +
  labs(x = "Null-model VPC (%)", y = "PCV (%)") +
  theme_bw(base_size = 12, base_family = "serif")
```
</details>

<details>
<summary><b>5. Tại sao PCV của scenario C trong Python là NaN còn trong R là 78%?</b></summary>

Trong Python, `σ²_null = 0` ở scenario C → `pcv(0, 0) = NaN`. Trong R, `σ²_null = 0.82` do khác biệt RNG → PCV tính được. **Cả hai đều hợp lệ.** Đây là ví dụ điển hình cho thấy kết quả VPC/PCV có thể dao động — và đó chính là điểm mà repo này muốn minh họa.
</details>

<details>
<summary><b>6. Tôi chỉ biết R, có cần Python không?</b></summary>

Không. R package `imaihda` là bản tái lập **hoàn chỉnh và độc lập**. Bạn có thể cài đặt, chạy toàn bộ pipeline, và sinh kết quả chỉ với R mà không cần Python.
</details>

---

## Tài liệu tham khảo

1. Evans CR, Williams DR, Onnela JP, Subramanian SV. A multilevel approach to modeling health inequalities at the intersection of multiple social identities. *SSM - Population Health*. 2018;6:149–157.
2. O'Sullivan JL, Alonso-Perez E, et al. Onset of Type 2 diabetes in adults aged 50 and older in Europe: an intersectional multilevel analysis of individual heterogeneity and discriminatory accuracy. *Diabetology & Metabolic Syndrome*. 2024;16:293. doi:10.1186/s13098-024-01533-3.
3. Elff M, Heisig JP, Schaeffer M, Shikano S. Multilevel analysis with few clusters. *British Journal of Political Science*. 2021;51:412–426.

---

## Giấy phép

MIT License — xem [LICENSE](LICENSE).

---

<a id="english"></a>
##  English

> [Tiếng Việt ở trên](#i-maihda-hic-mic-simulation-v31)

### What is this?

A synthetic-data stress-test workflow demonstrating that **VPC and PCV** — core I-MAIHDA metrics — are sensitive to outcome prevalence, sparse strata, and SES-patterned under-detection. A HIC-MIC difference in VPC/PCV does not necessarily indicate a different intersectional structure. Includes both **Python** (primary) and **R package `imaihda`** (full reproduction, installable).

### Quick start (R)

```r
remotes::install_github("nguyenminh2301/-i-maihda", subdir = "imaihda")
library(imaihda)

df <- simulate_intersectional_data(n = 2000, seed = 42)
res <- fit_imaihda(df)
res$vpc_null  # VPC from null model (%)
res$pcv       # PCV (%)
```

### Quick start (Python)

```bash
pip install -r requirements.txt
python scripts/run_all.py
# outputs: outputs/*.csv, figures/*.png
```

### Key finding

> A MIC cohort showing higher VPC or lower PCV than a HIC cohort does **not** necessarily mean the intersectional structure of inequality is different. The observed difference may reflect prevalence, sparse strata, or differential detection — not structural inequality.

### Cross-language validation

Both Python (PCG64) and R (Mersenne Twister) implementations pass **all 6/6 methodological benchmarks** despite different RNG engines. Numerical values differ, statistical conclusions match.

| Benchmark | Python | R |
|-----------|:------:|:--:|
| Additive-dominant detection | ✅ | ✅ |
| Interaction ↑ VPC | ✅ | ✅ |
| Interaction leaves residual | ✅ | ✅ |
| Detection ↓ prevalence | ✅ | ✅ |
| Detection masks VPC | ✅ | ✅ |
| Sparse strata flagged | ✅ | ✅ |

---

*Repository maintained by Minh Thien Nguyen. Last updated: June 2026.*
