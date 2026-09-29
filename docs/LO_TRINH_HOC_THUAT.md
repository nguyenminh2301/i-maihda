# Lộ trình nghiên cứu và phát triển dự án i-MAIHDA: từ minh chứng phương pháp đến so sánh thực nghiệm HIC–MIC có hiệu chỉnh sai lệch phát hiện

*Bản học thuật. Bản ngôn ngữ phổ thông có cùng nội dung: [LO_TRINH_DE_HIEU.md](LO_TRINH_DE_HIEU.md).*
*Cập nhật: 29/09/2026. Trạng thái tham chiếu: `imaihda` v0.4.0 (R) và `imaihda_sim` (Python), commit `d74ae01`.*

---

## Tóm tắt

**Bối cảnh.** Phân tích đa tầng giao thoa về tính dị biệt cá nhân và độ chính xác phân biệt (Intersectional Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy, I-MAIHDA) tóm tắt bất bình đẳng sức khoẻ giao thoa bằng hệ số phân chia phương sai (Variance Partition Coefficient, VPC) và tỉ lệ thay đổi phương sai (Proportional Change in Variance, PCV). Hai đại lượng này ngày càng được dùng để so sánh giữa quốc gia thu nhập cao (high-income countries, HIC) và thu nhập trung bình (middle-income countries, MIC).
**Vấn đề.** Các giai đoạn 1–3 của repo này cho thấy, trên dữ liệu mô phỏng, VPC và PCV biến thiên theo tỉ lệ hiện mắc, độ thưa của strata và sai lệch phát hiện có khuôn mẫu theo vị thế kinh tế – xã hội (SES-patterned differential detection), ngay cả khi cấu trúc giao thoa thật không đổi. Repo cũng đã xây dựng các công cụ phân tích sai lệch định lượng (quantitative bias analysis, QBA), một ước lượng hiệu chuẩn đồng thời và cận định danh từng phần sắc (sharp partial-identification bounds).
**Mục tiêu tối thượng.** Thiết lập một khung suy luận cho phép so sánh HIC–MIC về cấu trúc bất bình đẳng sức khoẻ giao thoa ở người trung–cao tuổi mà kết luận có hiệu lực dưới các giả định được nêu tường minh, và áp dụng khung đó lên dữ liệu cohort lão hoá hài hoà trong khuôn khổ luận án tiến sĩ.
**Lộ trình.** Sáu gói công việc (WP0–WP5) trong khoảng 18 tháng dẫn tới sáu sản phẩm: một gói phần mềm v1.0, ba bài báo (phương pháp–mô phỏng, lý thuyết định danh, ứng dụng thực nghiệm), một bảng kiểm báo cáo và luận án.

---

## 1. Bối cảnh và vấn đề nghiên cứu

### 1.1. Các đại lượng cốt lõi

Với cá thể $i$ thuộc stratum giao thoa $j \in \{1,\dots,K\}$ ($K = 36$ = giới 2 × học vấn 3 × tài sản 3 × nơi cư trú 2), mô hình logistic hiệu ứng ngẫu nhiên có dạng $\operatorname{logit} P(y^*_{ij}=1) = \beta_0 + \mathbf{x}_j^\top\boldsymbol\beta + u_j$, $u_j \sim N(0,\sigma^2_u)$. Trên thang biến tiềm ẩn:

$$\mathrm{VPC} = \frac{\sigma^2_u}{\sigma^2_u + \pi^2/3}, \qquad \mathrm{PCV} = \frac{\sigma^2_{\text{null}} - \sigma^2_{\text{main}}}{\sigma^2_{\text{null}}},$$

trong đó $\sigma^2_{\text{null}}$ là phương sai giữa strata của mô hình rỗng (chỉ có strata), còn $\sigma^2_{\text{main}}$ là phương sai còn lại sau khi đưa các hiệu ứng chính cộng tính vào mô hình (Evans và cs., 2018).

### 1.2. Vấn đề khả năng so sánh (comparability)

Kết quả quan sát là $y = y^* \cdot D$, với $D$ là chỉ báo phát hiện (chẩn đoán). Khi $P(D=1 \mid S=j)$ biến thiên theo SES, và mức biến thiên này khác nhau giữa các bối cảnh (thường lớn hơn ở MIC), thì hiệu $\mathrm{VPC}_{\text{MIC}} - \mathrm{VPC}_{\text{HIC}}$ tính trên $y$ lẫn lộn hai thành phần: khác biệt thật về cấu trúc giao thoa và khác biệt về quá trình đo lường. Độ thưa của strata và tỉ lệ hiện mắc thấp tạo thêm sai lệch do co rút (shrinkage) trong các bộ ước lượng nhanh.

### 1.3. Khoảng trống tri thức

Rà soát PubMed sơ bộ của dự án (ghi trong `PHASE3_CAUSAL_IDENTIFICATION.md` §1) không tìm thấy nghiên cứu I-MAIHDA nào xử lý sai số đo lường của biến kết cục. Gói CRAN `MAIHDA` (Bulut, 2026) mở rộng bề rộng mô hình (Bayes, trọng số khảo sát, dọc) nhưng không có công cụ cho sai lệch phát hiện. **Nhận định này cần được kiểm chứng bằng tổng quan hệ thống trước khi đưa vào bản thảo** (WP2, WP3).

---

## 2. Mục tiêu

### 2.1. Mục tiêu tối thượng

Xây dựng và kiểm chứng một khung suy luận, gồm lý thuyết định danh, bộ ước lượng, phần mềm và hướng dẫn báo cáo, để so sánh cấu trúc bất bình đẳng sức khoẻ giao thoa giữa HIC và MIC sao cho:

(i) các nguồn sai lệch không phải cấu trúc (tỉ lệ hiện mắc, độ thưa, sai lệch phát hiện) được nhận diện và hiệu chỉnh khi có thể;
(ii) khi không thể định danh điểm, kết luận được biểu đạt dưới dạng tập định danh (identified set) kèm ngưỡng lật kết luận (tipping point);
(iii) khung được áp dụng lên dữ liệu cohort lão hoá hài hoà để trả lời câu hỏi thực nghiệm của luận án.

*Căn cứ suy luận mục tiêu:* `docs/VERSION_3.1_QC_REPORT.md` (exposé PhD, thiết kế so sánh HIC–MIC có cấu trúc, HRS là đối chứng HIC ngoài châu Âu, lão hoá sinh học là mở rộng cơ chế có điều kiện sau audit hài hoà dấu ấn sinh học); `scripts/benchmark/benchmark_final.R` ("WZB Berlin PhD project"); `python/imaihda_sim/fit.py` ("The PhD empirical work should use full random-intercept MAIHDA models"). Exposé gốc không có trong repo; mục tiêu cần được đối chiếu với exposé.

### 2.2. Mục tiêu cụ thể

- **O1 (phần mềm).** Hoàn thiện `imaihda` thành gói R v1.0 đạt chuẩn CRAN, ngang bằng tính năng với bản Python, có kiểm thử chéo ngôn ngữ.
- **O2 (phương pháp–mô phỏng).** Định lượng độ nhạy của VPC/PCV và hiệu năng của các công cụ QBA theo một giao thức mô phỏng được đặc tả trước theo khung ADEMP (Aims, Data-generating mechanisms, Estimands, Methods, Performance measures).
- **O3 (lý thuyết).** Hoàn tất lý thuyết định danh cho các phiếm hàm phân chia phương sai: chứng minh hình thức Mệnh đề 4, mở rộng cận sắc sang PCV và sang thiết kế strata tổng quát (K bất kỳ).
- **O4 (thực nghiệm).** Ước lượng và so sánh VPC/PCV giữa các cohort HIC và MIC, trong đó xác suất phát hiện theo stratum được ước lượng trực tiếp từ mẫu phụ có dấu ấn sinh học.
- **O5 (chuyển giao).** Đề xuất bảng kiểm báo cáo I-MAIHDA khi nghi ngờ phát hiện không đồng đều.

---

## 3. Câu hỏi nghiên cứu và estimand

| Mã | Câu hỏi nghiên cứu | Estimand |
|---|---|---|
| RQ1 | VPC/PCV quan sát sai lệch thế nào so với giá trị quần thể dưới ba cơ chế (tỉ lệ hiện mắc, độ thưa, phát hiện)? | $\mathrm{VPC}(\sigma^2)$ với $\sigma^2 = \operatorname{Var}(\eta_j)$ của logit stratum thật (estimand quần thể, xem `PHASE3` §3.1) |
| RQ2 | Dưới giả định nào VPC được định danh điểm, định danh tập, hoặc không định danh? | Tập định danh $\Theta_I \subseteq [0,1]$ |
| RQ3 | Các công cụ QBA hoạt động ra sao khi giả định bị vi phạm? | Sai lệch, RMSE, độ phủ, tiêu chí không-kém-hơn (`METHODS_NOTE` §3) |
| RQ4 | Cận định danh có mở rộng được sang PCV và sang K tuỳ ý không? | $[\underline{\mathrm{PCV}}, \overline{\mathrm{PCV}}]$ |
| RQ5 | Sau khi tính đến sai lệch phát hiện, khác biệt HIC–MIC về VPC/PCV còn tồn tại không, và có vững không? | $\Delta = \mathrm{VPC}_{\text{MIC}} - \mathrm{VPC}_{\text{HIC}}$ và tập định danh của $\Delta$ |

**Estimand so sánh (điểm mới đề xuất cho RQ5).** Nếu cận của từng cohort $c$ là $[L_c, U_c]$ (từ `vpc_partial_bounds` với biên $\delta^{\max}_c$ riêng cho từng cohort) thì, vì tập định danh được xác định riêng cho từng cohort,

$$\Delta \in [\,L_{\text{MIC}} - U_{\text{HIC}},\; U_{\text{MIC}} - L_{\text{HIC}}\,].$$

Dấu của $\Delta$ được định danh khi khoảng này không chứa 0. Tổng quát hoá `detection_tipping_point()` cho phép tìm cặp $(\delta_{\text{HIC}}, \delta_{\text{MIC}})$ nhỏ nhất làm đảo dấu $\Delta$. Đây là phiên bản so sánh của lập luận kiểu giá trị E (E-value; VanderWeele & Ding, 2017).

---

## 4. Hiện trạng: những gì đã được chứng minh (dữ liệu mô phỏng)

| Thành phần | Bằng chứng trong repo | Kết quả chính | Trạng thái |
|---|---|---|---|
| Minh chứng độ nhạy VPC/PCV (5 kịch bản A–E) | `README.md`, `scripts/validation/`, `figures/scenario_vpc_pcv.png` | 6/6 tiêu chí đạt ở cả Python và R | Hoàn tất |
| Bộ ước lượng nhanh (method of moments, không trọng số) | `imaihda/R/diagnostics.R`, `imaihda/inst/benchmark/` | Nhanh hơn GLMM 92–144 lần; chênh VPC trung bình < 1 điểm % | Hoàn tất (R); Python chưa đồng bộ |
| Đối chứng với CRAN `MAIHDA` | `imaihda/tests/testthat/test-maihda-crossval.R` | Thành phần phương sai trùng khớp đến 1e-6; PCV đến 1e-4 | Hoàn tất |
| Hiệu chỉnh phát hiện `correct_detection_bias` | `detection.R` / `detection.py` | VPC quan sát 11.7%, thật 15.9%, hiệu chỉnh 15.7% | Hoàn tất (R + Python) |
| Hiệu chỉnh strata thưa `sparse_strata_vpc` | `uncertainty.R` / `uncertainty.py` | Sai lệch 4.2 → 0.3 điểm %; độ phủ CI 10/12 | Hoàn tất (R + Python) |
| Nghiên cứu độ vững khi vi phạm giả định | `docs/METHODS_NOTE_ROBUSTNESS.md` | Đổ vỡ khi dạng hàm score sai (sai lệch +9.73) và khi phát hiện phụ thuộc kết cục (`severity_weight` ≈ 0.8–0.9); không nhánh nào của `sparse_strata_vpc` bị đổ vỡ | Hoàn tất (chỉ Python) |
| Lý thuyết định danh (Mệnh đề 1–4) | `docs/PHASE3_CAUSAL_IDENTIFICATION.md` §2 | MĐ 1–3 có chứng minh; MĐ 4 mới được xác lập bằng thực nghiệm số | Một phần |
| Ước lượng hiệu chuẩn đồng thời `joint_calibrated_vpc` | `python/imaihda_sim/causal.py` | Vùng capstone: sai lệch −16.26 (naive), −6.72 (chỉ hiệu chỉnh phát hiện), +0.87 (joint); độ phủ 80–95% | Chỉ Python |
| Cận sắc `vpc_partial_bounds` | `causal.py`, `figures/phase3_causal_results.csv` | Chứa giá trị thật 45/45, kể cả khi định danh điểm bất khả | Chỉ Python |
| Kiểm thử tự động | `imaihda/tests/`, `python/tests/` | 79 lệnh `expect_*` (R); 33 test (Python) | Hoàn tất |

---

## 5. Phân tích khoảng trống

### 5.1. Nợ kỹ thuật (đã xác minh trong mã nguồn)

1. **Không nhất quán bộ ước lượng.** `python/imaihda_sim/fit.py` (`_between_stratum_variance`) vẫn dùng phương sai có trọng số độ chính xác, trong khi R (`imaihda/R/diagnostics.R`) dùng phương sai không trọng số từ v0.2.1. Chính `README.md` (mục Sparse-Strata) đã ghi nhận độ lệch dai dẳng của công thức có trọng số. `joint_calibrated_vpc` lại hiệu chuẩn trên thống kê có trọng số (`PHASE3` §3.2). Cần một định nghĩa chuẩn duy nhất cho thống kê naive và lý do chọn nó.
2. **Bảng kết quả lỗi thời.** Bảng kịch bản R trong `README.md` được tính bằng bộ ước lượng v0.2.0 (README tự ghi chú); bảng Python "unchanged since v3.1".
3. **Sai lệch tài liệu.** URL cài đặt và trường `URL`/`BugReports` trong `imaihda/DESCRIPTION` trỏ tới `nguyenminh2301/-i-maihda`, trong khi remote hiện tại là `nguyenminh2301/i-maihda`; README nhắc `benchmark2.R` (không tồn tại) và `source("benchmark_final.R")` ở thư mục gốc (thực tế ở `scripts/benchmark/`); bảng so sánh ghi "51 automated testthat assertions" trong khi hiện có 79.
4. **Tài liệu hàm R.** `imaihda/man/` nằm trong `.gitignore`, nên gói cài từ GitHub không có trang trợ giúp và `R CMD check` sẽ cảnh báo.
5. **Siêu dữ liệu tác giả.** `DESCRIPTION` liệt kê "Factory Droid" với vai trò `ctb`. Nên chuyển sang tuyên bố sử dụng công cụ AI theo chính sách của CRAN và của tạp chí.
6. **Tính tái lập.** `python/requirements.txt` không ghim phiên bản; R chưa có `renv.lock`.

### 5.2. Khoảng trống khoa học

1. **Ngang bằng R.** `joint_calibrated_vpc`, `vpc_partial_bounds` và các bộ sinh vi phạm giả định (`robustness.py`) chưa có bản R.
2. **Độ chính xác Monte Carlo.** Các nghiên cứu hiện dùng $n_{\text{rep}} = 12$–30. Với độ phủ danh nghĩa 95%, sai số chuẩn Monte Carlo (MCSE) của độ phủ là $\sqrt{0.95 \cdot 0.05 / n_{\text{rep}}} \approx 4.9$ điểm % khi $n_{\text{rep}} = 20$. Để MCSE ≤ 1 điểm % cần $n_{\text{rep}} \ge 475$. Với RMSE ≈ 9 điểm % ở vùng thưa, MCSE của sai lệch ≤ 0.5 điểm % đòi hỏi $n_{\text{rep}} \approx (9/0.5)^2 \approx 324$.
3. **Độ phủ dưới danh nghĩa ở vùng thưa** (80–85%), một phần do $\operatorname{Var}(\hat m)$ được coi là cố định trong CI của ước lượng joint (`PHASE3` §6).
4. **Mệnh đề 4** (tính đơn điệu toàn cục của ánh xạ hiệu chuẩn) mới là phỏng đoán có bằng chứng số.
5. **Cận cho PCV** và **tổng quát hoá K** chưa bắt đầu.
6. **Trọng số khảo sát và thiết kế dọc** chưa được hỗ trợ trong các công cụ QBA.
7. **Chưa có kiểm chứng ngoài** (external validation) trên dữ liệu thật.
8. **Trích dẫn chưa xác minh.** Khung định danh từng phần (Manski) và định danh dữ liệu thiếu bằng đồ thị (Mohan & Pearl) đang được nêu mà chưa có trích dẫn đã xác minh (`PHASE3` §8).

---

## 6. Kết quả kỳ vọng cuối cùng và tiêu chí nghiệm thu

| Mã | Sản phẩm | Tiêu chí nghiệm thu (Definition of Done) |
|---|---|---|
| **D1** | Gói `imaihda` v1.0.0 (R, CRAN) + `imaihda_sim` (Python) | `R CMD check --as-cran`: 0 lỗi, 0 cảnh báo; bảng ngang bằng tính năng R↔Python đầy đủ; các hàm tất định cho kết quả trùng khớp (dung sai 1e-8) trên cùng dữ liệu tổng hợp theo stratum; hàm ngẫu nhiên đồng nhất về phân phối; vignette "thang giả định"; bản lưu trữ có DOI |
| **D2** | Bài báo A: độ nhạy VPC/PCV và QBA (phương pháp + mô phỏng) | Giao thức ADEMP được đặc tả trước và lưu kèm mã; MCSE được báo cáo cho mọi chỉ số hiệu năng; $n_{\text{rep}}$ được biện minh; so sánh với GLMM đầy đủ; bản thảo đã nộp |
| **D3** | Bài báo B: định danh và định danh từng phần cho phiếm hàm phân chia phương sai | Chứng minh hình thức MĐ 1–3; MĐ 4 được chứng minh, hoặc được nêu tường minh là phỏng đoán kèm bằng chứng; cận PCV có kiểm chứng tối ưu toàn cục bằng liệt kê đầy đủ với K nhỏ; hỗ trợ K tuỳ ý; bản thảo đã nộp |
| **D4** | Bài báo C: so sánh HIC–MIC có hiệu chỉnh sai lệch phát hiện (chương thực nghiệm cốt lõi) | Kế hoạch phân tích được tiền đăng ký; báo cáo đủ các bậc của thang giả định; kiểm chứng ngoài của hiệu chỉnh bằng dấu ấn sinh học; tuân thủ STROBE; bản thảo đã nộp |
| **D5** | Bảng kiểm báo cáo I-MAIHDA dưới sai lệch phát hiện | Được áp dụng thử trên D4; công bố kèm D1 (vignette) hoặc D2 (phụ lục) |
| **D6** | Luận án | Ánh xạ chương: tổng quan → D2 → D3 → D4 → thảo luận tổng hợp và D5 |

---

## 7. Các gói công việc

### WP0. Củng cố nền tảng (tháng 1)

- **Nhiệm vụ:** xử lý toàn bộ mục 5.1; chọn và ghi tài liệu cho định nghĩa chuẩn của thống kê naive; chạy lại bảng kịch bản cho cả hai ngôn ngữ; commit trang `man/` (hoặc sinh tự động trong quy trình kiểm tra cục bộ); ghim phiên bản phụ thuộc (`requirements.txt` có phiên bản, `renv.lock`); tạo một lệnh kiểm tra cục bộ duy nhất (ví dụ `make check` chạy `pytest -q`, `devtools::test()` và `R CMD check --as-cran`).
- **Song song:** nộp hồ sơ xin quyền truy cập dữ liệu cohort (mục WP4), vì đây là đường găng (critical path) của toàn dự án.
- **Đầu ra:** v0.4.1; README (EN + VI) khớp mã; toàn bộ kiểm thử đạt trên máy cục bộ.

### WP1. Ngang bằng R và phát hành v1.0 (tháng 2–4)

- **Nhiệm vụ:** port `python/imaihda_sim/causal.py` (`joint_calibrated_vpc`, `vpc_partial_bounds` cùng các bộ tối ưu hoá `max` theo trao đổi toạ độ và `min` theo tìm kiếm tam phân) và `robustness.py` sang R; viết kiểm thử tương đương (bao gồm kiểm chứng tối ưu toàn cục bằng liệt kê $2^K$ với K = 3–11, như đã làm trong Python); bổ sung hỗ trợ trọng số khảo sát cho tỉ lệ hiện mắc theo stratum trong `correct_detection_bias` (cần cho WP4); viết vignette; gửi CRAN.
- **Tái sử dụng:** các hàm R hiện có `correct_detection_bias`, `vpc_detection_bounds`, `detection_tipping_point`, `sparse_strata_vpc`, `simulate_intersectional_data`, `between_stratum_variance`.
- **Đầu ra:** D1.

### WP2. Bài báo A: phương pháp và mô phỏng (tháng 4–7)

- **Aims:** RQ1, RQ3.
- **Data-generating mechanisms:** lưới nhân tố gồm $n \in \{2000, 3500, 12000\}$, dịch chuyển tỉ lệ hiện mắc $\{-2.1, -3, -4, -5\}$, `interaction_sd`, `detection_strength` $\delta \in [0, 1.2]$, `sparse`, cùng các nhánh vi phạm giả định M1-1 đến M2-3 đã có trong `python/imaihda_sim/robustness.py` và `scripts/validation/`.
- **Estimands:** VPC quần thể và PCV (tách rõ khỏi giá trị "refit trên `y_true`" hữu hạn mẫu, theo `PHASE3` §3.1).
- **Methods:** naive nhanh, GLMM đầy đủ (`glmer`), `correct_detection_bias`, `sparse_strata_vpc`, `joint_calibrated_vpc`, `vpc_partial_bounds`. Đề xuất thêm một đối chứng dựa trên mô hình: mô hình phân loại sai Bayes (misclassification model) có biến kết cục thật tiềm ẩn.
- **Performance measures:** sai lệch, RMSE, độ phủ, độ rộng CI, tỉ lệ cận chứa giá trị thật, mỗi chỉ số kèm MCSE; $n_{\text{rep}} \ge 500$ cho các kịch bản chính (mục 5.2).
- **Hạ tầng tính toán:** song song hoá trên các lõi CPU cục bộ; ghi hạt giống ngẫu nhiên theo từng replicate.
- **Đầu ra:** D2.

### WP3. Bài báo B: lý thuyết định danh (tháng 5–10, chồng lấn WP2)

- **MĐ 4, hướng chứng minh đề xuất.** Viết $\eta = m + \sigma Z$ với $Z$ đối xứng. Khi $\sigma_1 < \sigma_2$ thì $\sigma_1 Z \le_{cx} \sigma_2 Z$ (thứ tự lồi). Nếu $h(x) = E[T \mid \eta = x]$ lồi theo $x$ thì $E[T \mid \sigma^2]$ đơn điệu tăng. Việc cần chứng minh là tính lồi (hoặc một điều kiện yếu hơn) của $h$ dưới phép làm trơn Laplace và lấy mẫu nhị thức; nếu không được, dùng lập luận ghép cặp (coupling) theo từng stratum. Nếu vẫn không chứng minh được, giữ MĐ 4 là phỏng đoán có bằng chứng số (1.500 replicate/điểm lưới).
- **Cận PCV.** Với vector logit stratum $x$ trong hộp $[a,b]$: $\mathrm{PCV}(x) = 1 - \dfrac{x^\top (I-P)\,x}{x^\top (I-J)\,x}$, trong đó $P$ là phép chiếu lên không gian hiệu ứng chính và $J$ là phép chiếu lên hằng số. Đây là bài toán quy hoạch phân thức với các dạng toàn phương, không lồi. Hướng đề xuất: phương pháp tham số Dinkelbach (tìm nghiệm $\lambda$ của $\max_x [\,x^\top A x - \lambda\, x^\top B x\,] = 0$ bằng chia đôi), với bài toán con là quy hoạch toàn phương bất định trên hộp, giải toàn cục bằng rẽ nhánh–cận (branch-and-bound) cho K = 36, và dùng nới lỏng SDP để có chứng nhận cận. Cần báo cáo trung thực khả năng cận PCV trở nên không thông tin khi biên $\delta^{\max}$ rộng.
- **Tổng quát hoá K.** Viết lại bộ sinh dữ liệu cho thiết kế giai thừa tuỳ ý (số chiều và số mức bất kỳ), sau đó kiểm tra lại toàn bộ bộ kiểm thử.
- **Tổng quan tài liệu:** tổng quan hệ thống về xử lý sai số đo lường trong I-MAIHDA để kiểm chứng nhận định ở mục 1.3; xác minh các trích dẫn liệt kê ở mục 11.
- **Đầu ra:** D3.

### WP4. Ứng dụng thực nghiệm HIC–MIC (tháng 7–15)

- **Nguồn dữ liệu ứng viên** (cần đối chiếu với exposé): HRS (đã nêu trong báo cáo QC v3.1), SHARE, ELSA cho nhóm HIC; CHARLS, LASI, MHAS, ELSI-Brazil cho nhóm MIC. Ưu tiên các bộ biến đã được hài hoà sẵn (ví dụ Gateway to Global Aging Data) và kiểm tra lại từng biến.
- **Audit hài hoà (bước đi/không đi, go/no-go):** định nghĩa thống nhất cho giới, học vấn (ISCED rút gọn 3 mức), tài sản (tam phân vị trong từng quốc gia), nơi cư trú; kiểm tra phân bố cỡ strata để dự báo mức độ thưa; audit dấu ấn sinh học (loại mẫu: máu tĩnh mạch hay giọt máu khô; ngưỡng; hiệu chuẩn giữa các cohort), đúng như báo cáo QC v3.1 đã đặt điều kiện.
- **Kết cục ứng viên:** đái tháo đường (tự báo cáo được bác sĩ chẩn đoán so với HbA1c ≥ 6.5%) và tăng huyết áp (tự báo cáo so với huyết áp đo được hoặc đang dùng thuốc).
- **Ước lượng trực tiếp xác suất phát hiện.** Trong mẫu phụ có dấu ấn sinh học, $\hat d_j = \hat P(D=1 \mid y^*=1, S=j)$ (có trọng số khảo sát) chính là đại lượng $d_j$ của Mệnh đề 1. Giả định cần thiết là mẫu phụ đại diện cho stratum. Nhờ đó $\delta$ chuyển từ tham số độ nhạy thành tham số ước lượng được, với sai số lan truyền bằng bootstrap.
- **Thiết kế kiểm chứng ngoài.** So sánh (a) VPC tính trên kết cục định nghĩa bằng dấu ấn sinh học (không có sai lệch phát hiện nhưng mẫu nhỏ hơn và strata thưa hơn, nên dùng `joint_calibrated_vpc`/`sparse_strata_vpc`) với (b) VPC tính trên kết cục "được chẩn đoán" ở toàn mẫu, đã hiệu chỉnh bằng $\hat d_j$. Mức đồng thuận giữa (a) và (b) là bằng chứng thực nghiệm đầu tiên cho hiệu lực của công cụ. Đồng thời kiểm tra MĐ 2 bằng cách xem xác suất phát hiện có phụ thuộc mức độ nặng (ví dụ mức HbA1c trong nhóm bệnh) hay không.
- **Mô hình chính:** GLMM I-MAIHDA chuẩn bằng CRAN `MAIHDA` (engine `wemix` khi cần trọng số khảo sát); hiệu chỉnh tuổi (tuổi đưa vào cả mô hình rỗng và mô hình hiệu ứng chính, hoặc chuẩn hoá tuổi), vì cấu trúc tuổi khác nhau giữa HIC và MIC có thể gây nhiễu so sánh VPC.
- **Phân tích độ nhạy theo thang giả định** (`PHASE3` §5): bậc 1 là hiệu chỉnh điểm với $\hat d_j$; bậc 2 là hiệu chỉnh đồng thời khi strata thưa; bậc 3 là cận sắc và ngưỡng lật cho $\Delta$ (mục 3).
- **Quy trình:** tiền đăng ký trên OSF trước khi truy cập biến kết cục; báo cáo theo STROBE.
- **Quản trị dữ liệu:** toàn bộ dữ liệu có hạn chế được xử lý **cục bộ**, không tải lên dịch vụ đám mây hay công cụ AI trực tuyến; repo công khai chỉ chứa mã, dữ liệu mô phỏng và kết quả tổng hợp đạt ngưỡng công bố của từng thoả thuận dữ liệu.
- **Đầu ra:** D4.

### WP5. Tổng hợp và chuyển giao (tháng 15–18)

- **Nhiệm vụ:** hoàn thiện bảng kiểm D5 (các mục dự kiến: nêu cơ chế phát hiện giả định; báo cáo phân bố cỡ strata và số strata không có ca bệnh; báo cáo VPC/PCV thô cùng giá trị hiệu chỉnh; báo cáo cận và ngưỡng lật cho mọi so sánh giữa quần thể); hướng dẫn sử dụng; viết luận án.
- **Đầu ra:** D5, D6.

---

## 8. Lịch trình (mốc tháng 1 = 10/2026)

| Gói công việc | Q4/2026 | Q1/2027 | Q2/2027 | Q3/2027 | Q4/2027 | Q1/2028 | Sản phẩm |
|---|:-:|:-:|:-:|:-:|:-:|:-:|---|
| WP0 · Củng cố; nộp hồ sơ dữ liệu | ■ | | | | | | v0.4.1 |
| WP1 · Ngang bằng R, v1.0 | ■ | ■ | | | | | D1 |
| WP2 · Bài báo A | | ■ | ■ | | | | D2 |
| WP3 · Bài báo B | | ■ | ■ | ■ | | | D3 |
| WP4 · Ứng dụng HIC–MIC | | | ■ | ■ | ■ | | D4 |
| WP5 · Tổng hợp, luận án | | | | | ■ | ■ | D5, D6 |

**Mốc kiểm soát:** M1 (cuối Q4/2026): v0.4.1 và hồ sơ dữ liệu đã nộp. M2 (Q1/2027): D1 gửi CRAN. M3 (Q2/2027): D2 nộp; quyết định go/no-go của audit hài hoà. M4 (Q3/2027): D3 nộp. M5 (Q4/2027): D4 nộp. M6 (Q1/2028): bản thảo luận án.

---

## 9. Rủi ro và biện pháp giảm thiểu

| Rủi ro | Khả năng | Tác động | Biện pháp |
|---|:-:|:-:|---|
| Cấp quyền dữ liệu chậm | Cao | Cao | Nộp hồ sơ ở WP0; WP2–WP3 không phụ thuộc dữ liệu thật |
| Biến không hài hoà được giữa các cohort | Trung bình | Cao | Audit go/no-go sớm; thu hẹp tập cohort; dùng định nghĩa thô hơn nhưng nhất quán |
| Không có hoặc không đủ dữ liệu dấu ấn sinh học | Trung bình | Trung bình | Quay về $\delta$ làm tham số độ nhạy và cận sắc (công cụ đã có) |
| Mẫu phụ dấu ấn sinh học không đại diện cho stratum | Trung bình | Trung bình | Dùng trọng số mẫu phụ; phân tích độ nhạy với giả định MAR; báo cáo cận |
| Phát hiện phụ thuộc mức độ nặng (vi phạm MĐ 1) | Trung bình | Cao | Kiểm định trong mẫu phụ; nếu có, chuyển sang bậc 3 (cận) theo MĐ 2 |
| Không chứng minh được MĐ 4 | Trung bình | Thấp | Nêu là phỏng đoán có bằng chứng số; phép chiếu đẳng trương đã bảo đảm phép nghịch đảo xác định |
| Cận PCV không thông tin | Cao | Trung bình | Báo cáo trung thực; cung cấp cận có điều kiện theo $\delta^{\max}$ hẹp dần |
| Nhận định "khoảng trống tri thức" không đúng | Thấp | Trung bình | Tổng quan hệ thống ở WP3 trước khi viết phần giới thiệu |
| Chi phí tính toán mô phỏng lớn | Trung bình | Thấp | Bộ ước lượng nhanh; song song hoá cục bộ; ưu tiên kịch bản theo MCSE |

---

## 10. Nguyên tắc tái lập và quản trị dữ liệu

1. **Cục bộ trước hết.** Mọi phân tích trên dữ liệu có hạn chế chạy trên máy cục bộ. Không dùng dịch vụ đám mây, không đưa dữ liệu thật vào công cụ AI trực tuyến.
2. **Tách bạch nội dung repo.** Repo công khai chỉ gồm mã, dữ liệu mô phỏng, kết quả tổng hợp và tài liệu.
3. **Tái lập hoàn toàn.** Mỗi con số trong bản thảo được sinh ra bởi một script trong `scripts/`, với hạt giống ngẫu nhiên và phiên bản phụ thuộc được ghim.
4. **Minh bạch.** Giao thức mô phỏng (ADEMP) và kế hoạch phân tích thực nghiệm (OSF) được đặc tả trước. Mọi sai khác so với kế hoạch đều được báo cáo.
5. **Liêm chính khoa học.** Không có số liệu giả định trong tài liệu. Các nhận định chưa kiểm chứng được đánh dấu rõ ("cần xác minh").

---

## 11. Tài liệu tham khảo

**Đã có trong repo (đã trích dẫn trong README và `docs/`):**

1. Evans CR, Williams DR, Onnela J-P, Subramanian SV. A multilevel approach to modeling health inequalities at the intersection of multiple social identities. *SSM – Population Health*. 2018;6:149–157. doi:10.1016/j.ssmph.2018.08.005
2. O'Sullivan JL, Alonso-Perez E, et al. Onset of Type 2 diabetes in adults aged 50 and older in Europe: an intersectional multilevel analysis of individual heterogeneity and discriminatory accuracy. *Diabetology & Metabolic Syndrome*. 2024;16:293. doi:10.1186/s13098-024-01533-3
3. Elff M, Heisig JP, Schaeffer M, Shikano S. Multilevel analysis with few clusters: improving likelihood-based methods to provide unbiased estimates and accurate inference. *British Journal of Political Science*. 2021;51(1):412–426. doi:10.1017/S0007123419000097
4. Bulut O. MAIHDA: Intersectional Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy. R package version 0.1.11. 2026. https://cran.r-project.org/package=MAIHDA
5. VanderWeele TJ, Ding P. Sensitivity analysis in observational research: introducing the E-value. *Annals of Internal Medicine*. 2017;167(4):268–274. doi:10.7326/M16-2607

**Cần xác minh đầy đủ thông tin thư mục trước khi trích dẫn chính thức (WP3):**

- Morris, White & Crowther: khung ADEMP cho nghiên cứu mô phỏng (*Statistics in Medicine*, 2019).
- Manski: định danh từng phần (partial identification).
- Mohan & Pearl: định danh dữ liệu thiếu bằng mô hình đồ thị.
- Dinkelbach: quy hoạch phân thức phi tuyến (1967).
- Gateway to Global Aging Data: bộ dữ liệu hài hoà cho các cohort lão hoá.
- Tài liệu mô tả dữ liệu và dấu ấn sinh học của từng cohort ứng viên (HRS, SHARE, ELSA, CHARLS, LASI, MHAS, ELSI-Brazil).

---

*Tài liệu định hướng nghiên cứu. Mọi kết quả nêu ở mục 4 đến từ dữ liệu mô phỏng và các script trong repo. Tài liệu không đưa ra kết luận thực nghiệm nào về bất kỳ quần thể nào.*
