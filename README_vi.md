# I-MAIHDA HIC-MIC Simulation v3.2

> **English:** [README.md](README.md)

Quy trình mô phỏng dữ liệu tổng hợp kiểm định độ nhạy của **VPC** và **PCV** — hai chỉ số thống kê tóm tắt cốt lõi của I-MAIHDA — trước tác động của tỉ lệ hiện mắc, strata thưa, và sai số phát hiện có khuôn mẫu SES. Repository gồm bản **Python** (chính) và **R package `imaihda`** (tái lập hoàn chỉnh, có thể cài đặt).

⚠️ **Không dùng dữ liệu thật.** Repository này chỉ sử dụng dữ liệu tổng hợp. Không có tuyên bố thực nghiệm nào về bất kỳ quần thể nào. Đây là minh chứng phương pháp luận.

---

## Mục lục

1. [Câu hỏi nghiên cứu](#câu-hỏi-nghiên-cứu)
2. [Phương pháp](#phương-pháp)
3. [Kịch bản](#kịch-bản)
4. [Kết quả benchmark](#kết-quả-benchmark)
5. [Hình minh họa](#hình-minh-họa)
6. [R package `imaihda`](#r-package-imaihda)
7. [So sánh với CRAN `MAIHDA`](#so-sánh-với-cran-maihda)
8. [Đối chứng song ngữ](#đối-chứng-song-ngữ)
9. [R package so với script độc lập](#r-package-so-với-script-độc-lập)
10. [FAQs](#faqs)
11. [Tài liệu tham khảo](#tài-liệu-tham-khảo)

---

## Câu hỏi nghiên cứu

Nếu một cohort thu nhập trung bình (MIC) cho thấy **VPC cao hơn** hoặc **PCV thấp hơn** so với cohort thu nhập cao (HIC), liệu điều đó có nhất thiết nghĩa là cấu trúc bất bình đẳng giao thoa khác biệt?

**Trả lời: Không.** VPC và PCV có thể dao động theo tỉ lệ hiện mắc, strata giao thoa thưa, và sai số phát hiện có khuôn mẫu SES — ngay cả khi cấu trúc giao thoa thực sự không đổi. So sánh HIC‑MIC thô về VPC/PCV đòi hỏi chẩn đoán đi kèm.

---

## Phương pháp

Quy trình mô phỏng cá thể lồng trong **36 strata giao thoa** xác định bởi giới tính (2) × học vấn (3) × tài sản (3) × nông thôn/thiếu nguồn lực (2). Tính toán chẩn đoán I-MAIHDA nhanh qua empirical-stratum logit và mô hình logistic hiệu ứng chính, với bộ ước lượng **method-of-moments** trừ đi nhiễu nhị thức kỳ vọng khỏi phương sai trọng số quan sát được của phần dư cấp strata. Đây là chẩn đoán mô phỏng, không thay thế GLMM hiệu ứng ngẫu nhiên đầy đủ trong nghiên cứu thực nghiệm.

### Công thức

**VPC** — Hệ số phân vùng phương sai trên thang latent logistic:

$$VPC = \frac{\sigma^2_{\text{stratum}}}{\sigma^2_{\text{stratum}} + \pi^2/3} \times 100\%$$

**PCV** — Tỉ lệ thay đổi phương sai từ mô hình null (chỉ có strata) sang mô hình hiệu ứng chính cộng-gộp:

$$PCV = \frac{\sigma^2_{\text{null}} - \sigma^2_{\text{main}}}{\sigma^2_{\text{null}}} \times 100\%$$

Trong đó:
- $\sigma^2_{\text{null}}$ = phương sai giữa các strata từ mô hình null (chỉ có strata giao thoa)
- $\sigma^2_{\text{main}}$ = phương sai giữa strata còn lại sau hiệu ứng chính cộng-gộp của giới tính, học vấn, tài sản và nông thôn
- $\pi^2/3 \approx 3,\!29$ = phương sai mức cá thể của phân phối logistic chuẩn

### Quy trình sinh dữ liệu

1. **Phân bổ strata.** Cá thể được gán vào strata với xác suất đồng đều (6000 cá thể / 36 strata ≈ 167 mỗi strata), hoặc với trọng số gamma trong kịch bản thưa.
2. **Dự báo tuyến tính cộng-gộp.** $\eta = \beta_0 + \beta_1 \cdot \text{giới tính} + \beta_2 \cdot \text{học vấn} + \beta_3 \cdot \text{tài sản} + \beta_4 \cdot \text{nông thôn}$, với $\beta_0 = -2,\!10$ (tỉ lệ hiện mắc nền ~23%).
3. **Tương tác giao thoa dư (tuỳ chọn).** Hiệu ứng tương tác có cấu trúc được thêm ở cấp strata, sau đó trung tâm hoá để trực giao với intercept.
4. **Sai số phát hiện (tuỳ chọn).** Ca bệnh thật ít có khả năng được ghi nhận hơn ở strata thiệt thòi: $\text{logit}(P(\text{phát hiện})) = 2,\!0 - \delta \cdot \text{học vấn} - \delta \cdot \text{tài sản} - 0,\!4\delta \cdot \text{nông thôn}$.

---

## Kịch bản

| Kịch bản | Mô tả | Tham số chính |
|:--------:|-------|---------------|
| **A** | Gradient xã hội thuần cộng-gộp, phát hiện đồng đều | Mặc định |
| **B** | Tương tác giao thoa thực sự, phát hiện đồng đều | `interaction_sd = 0,90` |
| **C** | Cấu trúc cộng-gộp với sai số phát hiện theo khuôn mẫu SES | `detection_strength = 0,80` |
| **D** | Tương tác giao thoa + sai số phát hiện SES | `interaction_sd = 0,90`, `detection_strength = 0,80` |
| **E** | Tương tác giao thoa, bệnh hiếm, strata thưa | `n = 3500`, `prevalence_shift = -3,00`, `interaction_sd = 0,90`, `sparse = TRUE` |

---

## Kết quả benchmark

### Ước lượng theo kịch bản

#### Python (PCG64 RNG, NumPy `default_rng`)

| Kịch bản | Tỉ lệ hiện mắc | VPC null | VPC main | PCV | Cỡ strata nhỏ nhất |
|:--------:|:--------------:|:--------:|:--------:|:---:|:------------------:|
| **A** | 23,3% | 4,32 | 0,00 | **100,0** | 144 |
| **B** | 27,1% | 22,58 | 15,78 | **35,8** | 144 |
| **C** | 11,3% | 0,00 | 0,00 | NaN | 144 |
| **D** | 13,7% | 13,68 | 8,80 | **39,1** | 144 |
| **E** | 9,1% | 14,70 | 9,44 | **39,5** | 1 |

#### R (package `imaihda`, Mersenne Twister RNG)

| Kịch bản | Tỉ lệ hiện mắc | VPC null | VPC main | PCV | Cỡ strata nhỏ nhất |
|:--------:|:--------------:|:--------:|:--------:|:---:|:------------------:|
| **A** | 23,6% | 4,50 | 0,00 | **100,0** | 130 |
| **B** | 26,3% | 17,20 | 9,12 | **51,7** | 130 |
| **C** | 11,6% | 0,82 | 0,18 | **78,0** | 130 |
| **D** | 13,0% | 16,06 | 8,02 | **54,5** | 130 |
| **E** | 11,6% | 19,47 | 2,98 | **87,3** | 2 |

### Tiêu chí pass/fail (cả hai ngôn ngữ cho kết quả giống hệt)

| # | Tiêu chí | Python | R |
|---|----------|:------:|:--:|
| 1 | A thuần cộng-gộp: PCV ≥ 80, VPC_main < 1 | ✅ | ✅ |
| 2 | B tương tác làm tăng VPC: VPC_null(B) > VPC_null(A) + 5 điểm phần trăm | ✅ | ✅ |
| 3 | B để lại phương sai dư: PCV < 70 | ✅ | ✅ |
| 4 | C sai số phát hiện làm giảm tỉ lệ hiện mắc quan sát | ✅ | ✅ |
| 5 | D sai số phát hiện che lấp VPC tương tác: VPC_null(D) < VPC_null(B) | ✅ | ✅ |
| 6 | E strata thưa được gắn cờ: min_n(E) < min_n(B) | ✅ | ✅ |

> **Kết luận:** Cả Python và R đều xác nhận rằng VPC và PCV dao động theo tỉ lệ hiện mắc, strata thưa, và sai số phát hiện. So sánh HIC‑MIC thô không thể diễn giải được nếu thiếu chẩn đoán strata đi kèm.

---

## Hình minh họa

### 1. Bản đồ VPC-PCV theo kịch bản

![VPC-PCV scenario plot](figures/scenario_vpc_pcv.png)

**Diễn giải:** Kịch bản A (góc trên bên trái) thể hiện cấu trúc thuần cộng-gộp (PCV = 100%). Thêm tương tác giao thoa thực sự (B) đẩy điểm sang phải (VPC cao hơn) và xuống dưới (PCV thấp hơn). Kịch bản D cho thấy sai số phát hiện có thể che lấp VPC ngay cả khi cùng một mức tương tác dư. Kịch bản E minh họa tác động của strata thưa lên cả VPC và PCV.

### 2. Quét sai số phát hiện

![Detection sweep](figures/detection_sweep.png)

**Diễn giải:** Khi cường độ sai số phát hiện theo khuôn mẫu SES tăng, tỉ lệ hiện mắc quan sát giảm đơn điệu (đường đứt nét). VPC thể hiện phản ứng **không đơn điệu**: ban đầu giảm (che lấp) và sau đó có thể tăng trở lại ở mức sai số cực cao — vì một số strata mất gần như toàn bộ ca quan sát trong khi các strata khác vẫn giữ được ca phát hiện. Tính không đơn điệu này nhấn mạnh lý do không thể bỏ qua sai số phát hiện trong so sánh VPC xuyên cohort.

---

## R package `imaihda`

Package R có thể cài đặt, có tài liệu đầy đủ, chứa toàn bộ quy trình mô phỏng và chẩn đoán. 6 hàm được export; 6 hàm nội bộ.

### Cài đặt

```r
# Từ GitHub
remotes::install_github("nguyenminh2301/-i-maihda", subdir = "imaihda")

# Hoặc clone về cài cục bộ
# git clone https://github.com/nguyenminh2301/-i-maihda.git
# devtools::install("đường-dẫn/-i-maihda/imaihda")
```

**Yêu cầu:** R ≥ 4.0. Phụ thuộc: `stats` (base R). Gợi ý: `ggplot2`, `testthat`, `viridis`.

### Sử dụng

```r
library(imaihda)
```

#### `vpc_latent()` — Tính VPC từ phương sai strata

```r
vpc_latent(0,5)       # 13,2% — bất bình đẳng giữa strata ở mức trung bình
vpc_latent(0)         # 0%
vpc_latent(pi^2 / 3)  # 50% — phương sai strata bằng phương sai cá thể
```

#### `pcv()` — Tính tỉ lệ thay đổi phương sai

```r
pcv(1,0, 0,25)  # 75% — phần lớn phương sai được giải thích bởi hiệu ứng cộng-gộp
pcv(0,5, 0,4)   # 20% — còn nhiều tương tác dư
pcv(0, 0)       # NaN — không xác định khi phương sai null ≤ 0
```

#### `simulate_intersectional_data()` — Sinh dữ liệu tổng hợp

```r
# Cơ bản (thuần cộng-gộp, phát hiện đồng đều)
df <- simulate_intersectional_data(n = 2000, seed = 42)

# Có tương tác giao thoa
df_b <- simulate_intersectional_data(n = 2000, interaction_sd = 0,9, seed = 42)

# Có sai số phát hiện theo SES
df_c <- simulate_intersectional_data(n = 2000, detection_strength = 0,8, seed = 42)

# Strata thưa, bệnh hiếm
df_e <- simulate_intersectional_data(
  n = 1000, prevalence_shift = -3,0,
  interaction_sd = 0,9, sparse = TRUE, seed = 42
)

# So sánh tỉ lệ hiện mắc quan sát và thực khi có sai số phát hiện
mean(df_c$y)       # quan sát (thấp hơn do phát hiện thiếu)
mean(df_c$y_true)  # thực (cao hơn)
```

#### `fit_imaihda()` — Chẩn đoán MAIHDA một lần gọi

```r
df  <- simulate_intersectional_data(n = 3000, seed = 123)
res <- fit_imaihda(df)

res$n_strata              # 36
res$overall_prevalence    # ~0,23
res$vpc_null              # VPC từ mô hình null (%)
res$vpc_main              # VPC sau hiệu ứng chính cộng-gộp (%)
res$pcv                   # Tỉ lệ thay đổi phương sai (%)
res$var_null              # Phương sai giữa strata (null)
res$var_main              # Phương sai giữa strata (main)
res$min_stratum_n         # Cỡ strata nhỏ nhất
```

#### `scenario_grid()` + `evaluate_benchmarks()` — Quy trình đầy đủ

```r
grid <- scenario_grid()
results <- do.call(rbind, lapply(names(grid), function(nm) {
  as.data.frame(fit_scenario(nm, grid[[nm]]))
}))
benchmarks <- evaluate_benchmarks(results)
print(benchmarks)  # 6 dòng pass/fail
```

#### Chạy kiểm định

```r
devtools::test("imaihda")   # 39 testthat assertions
```

---

## Đối chứng song ngữ

### Khác biệt về RNG

| Khía cạnh | Python | R |
|-----------|--------|---|
| **Engine** | PCG64 (`numpy.random.default_rng`) | Mersenne Twister (`set.seed`) |
| **Hạt giống** | 42 | 42 |
| **Kết quả số** | Khác | Khác |
| **Kết quả benchmark** | 6/6 pass | 6/6 pass |

### So sánh từng chỉ số

| Chỉ số | Python (điển hình) | R (điển hình) | Nhất quán |
|--------|:------------------:|:-------------:|:---------:|
| VPC_null(A) | 4,32 | 4,50 | ✅ Thấp ở cả hai |
| VPC_null(B) > VPC_null(A) | Có (22,58 > 4,32) | Có (17,20 > 4,50) | ✅ |
| PCV(A) | 100,0 | 100,0 | ✅ Thuần cộng-gộp |
| PCV(B) < 70 | Có (35,8) | Có (51,7) | ✅ Còn tương tác dư |
| Tỉ lệ hiện mắc C/A | 11,3/23,3 | 11,6/23,6 | ✅ Giảm ~50% |
| VPC(D) < VPC(B) | Có (13,68 < 22,58) | Có (16,06 < 17,20) | ✅ Hiệu ứng che lấp |
| Strata thưa ở E | min_n = 1 | min_n = 2 | ✅ Được gắn cờ |

> Cả hai bản triển khai đều đạt **kết luận định tính giống hệt nhau**. Khác biệt số liệu phát sinh từ khác biệt engine RNG và là điều được kỳ vọng trong bất kỳ tái lập song ngữ nào sử dụng mô phỏng ngẫu nhiên. Chúng không ảnh hưởng đến diễn giải khoa học.

---

## So sánh với CRAN `MAIHDA`

Package [`MAIHDA`](https://cran.r-project.org/package=MAIHDA) trên CRAN (Bulut 2026, v0.1.11) là **công cụ thực nghiệm chuẩn vàng** cho MAIHDA giao thoa. Package dùng ước lượng GLMM đầy đủ qua `lme4::glmer()`/`lmer()` hoặc `brms`, hỗ trợ trọng số khảo sát, bootstrap CI, và bảng điều khiển Shiny. Package `imaihda` này là **bộ công cụ mô phỏng và stress-test bổ trợ**.

### So sánh tính năng

| Khả năng | CRAN `MAIHDA` | `imaihda` |
|----------|:------------:|:---------:|
| VPC & PCV dựa trên GLMM | ✅ `lme4`/`brms` | ✅ `method="glmer"` |
| Chẩn đoán xấp xỉ nhanh | — | ✅ `method="fast"` |
| Mô phỏng dữ liệu tổng hợp | — | ✅ `simulate_intersectional_data()` |
| Sai số phát hiện theo khuôn mẫu SES | — | ✅ cường độ tùy chỉnh |
| Kịch bản stress-test dựng sẵn | — | ✅ 5 kịch bản (A–E) |
| Đánh giá benchmark tự động | — | ✅ 6 tiêu chí pass/fail |
| Song ngữ Python + R | — | ✅ hai bản triển khai |
| Khoảng tin cậy bootstrap | ✅ tham số | — |
| Trọng số khảo sát (design-weighted) | ✅ `WeMix` | — |
| Bảng điều khiển Shiny | ✅ `run_maihda_app()` | — |
| Phân rã PCV từng bước | ✅ `stepwise_pcv()` | — |
| So sánh nhóm | ✅ `compare_maihda_groups()` | — |
| Mô hình chéo/phân loại / dọc | ✅ | — |
| Độ chính xác phân biệt (AUC, MOR) | ✅ | — |

### Đối chứng chéo

Chúng tôi đã kiểm định `imaihda(method="glmer")` với CRAN `MAIHDA` trên dữ liệu NHANES đi kèm (`maihda_health_data`). Cả hai package cho ra **thành phần phương sai và ước lượng VPC/PCV giống hệt nhau**:

| Chỉ số | CRAN `MAIHDA` | `imaihda` (glmer) | Khớp |
|--------|:------------:|:-----------------:|:----:|
| Phương sai giữa strata (null) | 2,831 | 2,831 | ✅ 1e-6 |
| Phương sai giữa strata (main) | 0,492 | 0,492 | ✅ 1e-6 |
| VPC (null) | 0,0636 | 0,0636 | ✅ 1e-6 |
| PCV | 0,826 | 0,826 | ✅ 1e-4 |

> **Điểm mấu chốt:** CRAN `MAIHDA` báo cáo VPC từ mô hình REML (mặc định cho mô hình tuyến tính hỗn hợp) và PCV từ mô hình được refit bằng ML (qua `maihda_pcv_refit_ml()`). Với biến nhị phân (logistic GLMM), sự phân biệt này không liên quan vì `glmer()` luôn dùng ML. `method="glmer"` của chúng tôi dùng ML xuyên suốt cho biến nhị phân, khớp chính xác với ước lượng logistic của CRAN `MAIHDA`.

49 assertions testthat (gồm 12 test đối chứng chéo) xác nhận sự tương đương về số học. Xem `tests/testthat/test-maihda-crossval.R`.

### Khi nào dùng package nào

| Trường hợp sử dụng | Package khuyến nghị |
|--------------------|---------------------|
| MAIHDA thực nghiệm với dữ liệu khảo sát thật | CRAN `MAIHDA` |
| Stress-test phương pháp luận về độ nhạy VPC/PCV | `imaihda` |
| Sinh dữ liệu tổng hợp có sai số phát hiện | `imaihda` |
| Kiểm tra khả năng tái lập song ngữ (Python ↔ R) | `imaihda` |
| MAIHDA chuẩn publication với bootstrap CI | CRAN `MAIHDA` |
| Khám phá tương tác (Shiny) | CRAN `MAIHDA` |

## R package so với script độc lập

R package `imaihda` (v3.2) thay thế các script R độc lập trước đây (`R/*.R`, v3.1).

| Tiêu chí | Script độc lập (v3.1) | R package (v3.2) |
|----------|------------------------|-------------------|
| **Cấu trúc** | File `.R` rời, `source()` thủ công | Package chuẩn: DESCRIPTION, NAMESPACE |
| **Cài đặt** | Sao chép file, `source()` thủ công | `install_github()` hoặc `devtools::install()` |
| **Tài liệu** | Chỉ có comment nội bộ | Roxygen2 với `@examples`, `@references`, `@export` |
| **API xuất ra** | Không phân biệt public/private | 6 hàm xuất, 6 hàm nội bộ |
| **Kiểm định** | 4 khối `test_that` tạm | 39 assertions `testthat` tự động |
| **Tính khả chuyển** | Gắn với thư mục dự án WZB | Độc lập, dùng được ở mọi dự án |
| **Khả năng tái lập** | Cùng thuật toán | Cùng thuật toán — **kết quả giống 100%** ở cùng seed |

> **Tính nhất quán đã được xác nhận:** Package sử dụng cùng logic tính toán với script độc lập. Ở cùng hạt giống, kết quả số giống hệt từng bit vì thuật toán và lời gọi RNG không thay đổi — chỉ khác về cách tổ chức mã nguồn.

---

## FAQs

<details>
<summary><strong>1. Đây có phải là bộ ước lượng mới cho MAIHDA không?</strong></summary>

Không. Đây là **minh chứng phương pháp luận** sử dụng chẩn đoán empirical-logit nhanh để stress-test lặp lại. Nó không thay thế ước lượng GLMM hiệu ứng ngẫu nhiên đầy đủ (ví dụ: `lme4::glmer` trong R hoặc các bản triển khai mixed-model tương đương). Với nghiên cứu thực nghiệm, kết quả cần được kiểm chứng dựa trên chiến lược mô hình hoá của nhóm nghiên cứu mục tiêu.
</details>

<details>
<summary><strong>2. Tại sao Python và R cho ra số liệu khác nhau?</strong></summary>

Vì chúng dùng **bộ sinh số ngẫu nhiên khác nhau**: PCG64 trong NumPy và Mersenne Twister trong R. Cùng hạt giống (`42`) nhưng chuỗi sinh ra khác nhau, dẫn đến tập dữ liệu mô phỏng khác nhau và do đó ước lượng điểm VPC/PCV khác nhau. **Cả 6/6 benchmark đều pass ở cả hai ngôn ngữ**, và mọi kết luận định tính đều giống hệt. Đây là hành vi được kỳ vọng trong bất kỳ tái lập ngẫu nhiên song ngữ nào.
</details>

<details>
<summary><strong>3. Tôi có thể dùng package này với dữ liệu thật không?</strong></summary>

Có thể dùng `fit_imaihda()` để chẩn đoán thăm dò nhanh, nhưng hàm này sử dụng **xấp xỉ method-of-moments** trừ đi nhiễu nhị thức ước lượng — không phải là bộ ước lượng GLMM đầy đủ. Với công bố thực nghiệm, hãy dùng mô hình logistic hiệu ứng ngẫu nhiên phù hợp như `lme4::glmer(y ~ (1 | stratum) + covariates, family = binomial)`.
</details>

<details>
<summary><strong>4. Làm sao để vẽ lại các hình?</strong></summary>

```r
library(imaihda)
library(ggplot2)

grid <- scenario_grid()
results <- do.call(rbind, lapply(names(grid), function(nm) {
  as.data.frame(fit_scenario(nm, grid[[nm]]))
}))

ggplot(results, aes(vpc_null, pcv, label = scenario)) +
  geom_point(size = 3, color = "#21918c") +
  geom_text(hjust = -0,3, family = "serif") +
  labs(x = "VPC mô hình null (%)", y = "PCV (%)") +
  theme_bw(base_size = 12, base_family = "serif")
```
</details>

<details>
<summary><strong>5. Tại sao PCV ở kịch bản C của Python là NaN còn của R là 78%?</strong></summary>

Trong lần chạy Python, kịch bản C cho ra $\sigma^2_{\text{null}} = 0$, do đó `pcv(0, 0) = NaN` theo định nghĩa. Trong lần chạy R, $\sigma^2_{\text{null}} = 0,\!82$ do khác biệt RNG, nên PCV tính được. **Cả hai đều là kết quả hợp lệ.** Chính sự khác biệt này minh họa cho luận điểm của repository: ước lượng VPC/PCV có thể dao động giữa các lần thực hiện ngẫu nhiên, và các giá trị phương sai giữa strata nhỏ cần được diễn giải thận trọng.
</details>

<details>
<summary><strong>6. Tôi có cần Python không nếu chỉ làm việc với R?</strong></summary>

Không. R package `imaihda` là bản tái lập **hoàn chỉnh và độc lập**. Bạn có thể cài đặt, chạy toàn bộ quy trình, và tạo mọi kết quả chỉ với R.
</details>

---

## Tài liệu tham khảo

1. Evans CR, Williams DR, Onnela J-P, Subramanian SV. A multilevel approach to modeling health inequalities at the intersection of multiple social identities. *SSM - Population Health*. 2018;6:149–157. doi:10.1016/j.ssmph.2018.08.005
2. O'Sullivan JL, Alonso-Perez E, et al. Onset of Type 2 diabetes in adults aged 50 and older in Europe: an intersectional multilevel analysis of individual heterogeneity and discriminatory accuracy. *Diabetology & Metabolic Syndrome*. 2024;16:293. doi:10.1186/s13098-024-01533-3
3. Elff M, Heisig JP, Schaeffer M, Shikano S. Multilevel analysis with few clusters: improving likelihood-based methods to provide unbiased estimates and accurate inference. *British Journal of Political Science*. 2021;51(1):412–426. doi:10.1017/S0007123419000097

---

## Giấy phép

MIT — xem [LICENSE](LICENSE).

---

*Bảo trì bởi Minh Thien Nguyen. Cập nhật lần cuối: tháng 6 năm 2026.*
