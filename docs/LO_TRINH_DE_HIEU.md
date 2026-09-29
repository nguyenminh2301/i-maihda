# Lộ trình dự án i-MAIHDA — Bản dễ hiểu

> Bản này viết cho người không chuyên thống kê: đồng nghiệp, hội đồng, người hướng dẫn ở lĩnh vực khác, hoặc chính chúng ta khi cần nhìn lại toàn cảnh. Đọc khoảng 10 phút.
> Bản học thuật có cùng nội dung nhưng chặt chẽ hơn: [LO_TRINH_HOC_THUAT.md](LO_TRINH_HOC_THUAT.md).
>
> Cập nhật: 29/09/2026 · Dựa trên trạng thái repo ở phiên bản `imaihda` v0.4.0.

---

## 1. Dự án này để làm gì? (nói trong một câu)

**Giúp các nhà nghiên cứu so sánh công bằng mức độ bất bình đẳng sức khoẻ giữa nước giàu và nước thu nhập trung bình, khi hai nơi "đo bệnh" không giống nhau.**

### Một ví dụ đời thường

Hãy tưởng tượng ta muốn biết ở nước nào, *việc bạn thuộc nhóm xã hội nào* (nam hay nữ, học nhiều hay ít, giàu hay nghèo, ở thành phố hay nông thôn) ảnh hưởng mạnh hơn tới nguy cơ mắc bệnh, ví dụ tiểu đường ở người trên 50 tuổi.

Ta chia dân số thành **36 nhóm** (2 giới × 3 mức học vấn × 3 mức tài sản × 2 nơi ở) rồi xem tỉ lệ bệnh chênh nhau giữa các nhóm đến mức nào. Phương pháp chuẩn để làm việc này tên là **I-MAIHDA**, và nó cho ra hai con số chính:

- **VPC**: trong toàn bộ khác biệt về nguy cơ bệnh giữa từng người, bao nhiêu phần trăm nằm ở chuyện *người đó thuộc nhóm nào*. VPC càng cao thì "nhóm xã hội" càng quyết định sức khoẻ.
- **PCV**: trong phần khác biệt giữa các nhóm đó, bao nhiêu phần trăm giải thích được bằng cách *cộng dồn* từng yếu tố riêng (giới + học vấn + tài sản + nơi ở). Phần còn lại là "hiệu ứng giao thoa": những bất lợi chồng lên nhau theo kiểu không cộng đơn giản.

### Vấn đề: chiếc cân bị lệch

So sánh VPC của nước giàu với nước thu nhập trung bình giống như **cân hai giỏ trái cây bằng hai chiếc cân khác nhau**. Nếu một chiếc cân bị lệch, ta sẽ kết luận sai giỏ nào nặng hơn.

Dự án chỉ ra ba cách chiếc cân có thể bị lệch:

| Cái bẫy | Nói dễ hiểu | Hậu quả |
|---|---|---|
| **Tỉ lệ bệnh khác nhau** | Bệnh hiếm ở nơi này, phổ biến ở nơi kia | VPC và PCV thay đổi dù cấu trúc bất bình đẳng thật không đổi |
| **Nhóm quá ít người** ("strata thưa") | Có nhóm chỉ vài người, số liệu chập chờn | Công cụ tính nhanh đánh giá thấp mức chênh lệch thật |
| **Bệnh không được phát hiện đều** | Người nghèo, ít học, ở nông thôn ít được khám nên bệnh bị "giấu" | Chênh lệch nhìn thấy khác chênh lệch thật, có thể bị che đi hoặc phóng đại |

Cái bẫy thứ ba đáng lo nhất ở các nước thu nhập trung bình, vì ở đó việc được chẩn đoán phụ thuộc nhiều vào hoàn cảnh kinh tế.

**Câu trả lời chính của dự án đến nay:** VPC cao hơn hay PCV thấp hơn ở một nước *không tự động* có nghĩa là bất bình đẳng ở nước đó có cấu trúc khác. Phải kiểm tra ba cái bẫy trên trước.

---

## 2. Mục tiêu tối thượng (đích đến cuối cùng)

Dự án có đích đến rõ ràng hơn một bộ phần mềm. Theo các tài liệu trong repo (báo cáo QC phiên bản 3.1, ghi chú trong mã về "PhD empirical work", dự án PhD tại WZB Berlin), mục tiêu cuối là:

> **Tạo ra cách so sánh bất bình đẳng sức khoẻ giao thoa giữa nước thu nhập cao và thu nhập trung bình ở người trung–cao tuổi sao cho kết luận đáng tin: tách được khác biệt *thật* khỏi khác biệt *giả* do cách đo; rồi dùng cách đó trên dữ liệu thật của các cohort lão hoá để trả lời câu hỏi của luận án.**

Nói gọn, mục tiêu có ba lớp:

1. **Chứng minh vấn đề:** cho thấy so sánh "thô" có thể sai và sai theo cách nào. *(Đã làm xong phần lớn.)*
2. **Làm ra công cụ sửa:** công cụ vừa sửa được sai lệch vừa nói thật khi không sửa được. *(Đã có, đang hoàn thiện.)*
3. **Dùng công cụ trên dữ liệu thật:** trả lời câu hỏi thực nghiệm của luận án. *(Chưa bắt đầu.)*

> Lưu ý: mục tiêu trên được suy ra từ các tài liệu trong repo. Bản exposé gốc của luận án không nằm trong repo, nên cần đối chiếu lại với exposé.

---

## 3. Đến nay đã làm được gì?

Dự án đã qua ba giai đoạn. Tất cả kết quả dưới đây đến từ **dữ liệu mô phỏng** (dữ liệu giả lập, biết trước đáp án đúng), chưa dùng dữ liệu thật.

| Việc đã xong | Kết quả bằng lời thường |
|---|---|
| Chứng minh ba cái bẫy | 5 kịch bản mô phỏng (A–E), 6/6 phép kiểm tra đều đạt, ở cả Python và R |
| Công cụ tính nhanh | Nhanh hơn phương pháp chuẩn khoảng 92–144 lần; VPC chênh trung bình dưới 1 điểm phần trăm |
| Khớp với gói chuẩn quốc tế | Chế độ GLMM cho kết quả trùng với gói CRAN `MAIHDA` (sai khác không quá 0.0001) |
| Sửa lỗi "bệnh bị giấu" | Ví dụ: VPC nhìn thấy 11.7%, VPC thật 15.9%, sau khi sửa được 15.7% |
| Sửa lỗi "nhóm quá ít người" | Sai lệch trung bình giảm từ 4.2 xuống 0.3 điểm phần trăm |
| Thử công cụ trong điều kiện xấu | Biết rõ khi nào công cụ còn giúp ích, khi nào bắt đầu làm hại (xem mục 5) |
| Công cụ sửa cả hai lỗi cùng lúc | Sai lệch giảm từ −16.3 xuống +0.9 điểm phần trăm trong tình huống khó nhất |
| "Khoảng an toàn" khi không sửa được | Khoảng giá trị luôn chứa đáp án đúng: 45/45 lần thử, kể cả khi mọi cách sửa khác đều thất bại |
| Kiểm thử tự động | 79 phép kiểm trong R, 33 bài kiểm trong Python |

Chi tiết: [README_vi.md](../README_vi.md), [METHODS_NOTE_ROBUSTNESS.md](METHODS_NOTE_ROBUSTNESS.md), [PHASE3_CAUSAL_IDENTIFICATION.md](PHASE3_CAUSAL_IDENTIFICATION.md).

### Bài học quan trọng nhất đến nay

1. **Công cụ sửa "bệnh bị giấu" chịu được việc chọn thiếu một yếu tố**, nhưng **hỏng nếu giả định sai hình dạng** của quy luật phát hiện: khi đó nó sửa quá tay sang phía ngược lại.
2. **Nếu chính việc có bệnh làm người ta ít được phát hiện hơn** (không chỉ do nghèo), thì *không có cách nào* tính ra một con số đúng duy nhất. Khi đó cách trung thực là báo cáo **một khoảng** chứa giá trị thật.
3. **Không được áp hai cách sửa nối tiếp nhau.** Sửa lỗi này rồi sửa tiếp lỗi kia sẽ sửa quá tay (+8.8 điểm phần trăm). Phải dùng công cụ sửa đồng thời.

---

## 4. Còn thiếu gì?

| Hạng mục còn thiếu | Vì sao quan trọng |
|---|---|
| Một số lỗi nhỏ và chỗ không khớp trong tài liệu, mã | Người ngoài cài đặt hoặc đọc sẽ gặp lỗi, mất tin tưởng |
| Bản Python và bản R tính "cách nhanh" chưa giống nhau hoàn toàn | Hai ngôn ngữ phải cho cùng một kết luận khi dùng cùng công thức |
| Công cụ sửa đồng thời và "khoảng an toàn" mới có trong Python | Phần lớn người làm dịch tễ dùng R |
| Số lần lặp mô phỏng còn ít (20–30 lần) | Chưa đủ chính xác để đăng tạp chí phương pháp hàng đầu |
| Một số chứng minh toán học mới dừng ở kiểm tra bằng máy | Bài báo lý thuyết cần chứng minh chặt |
| "Khoảng an toàn" mới có cho VPC, chưa có cho PCV | PCV là con số quan trọng thứ hai của I-MAIHDA |
| Mới thử với đúng 36 nhóm | Nghiên cứu thật có thể có số nhóm khác |
| **Chưa dùng dữ liệu thật** | Đây là mục tiêu cuối của luận án |

---

## 5. Kết quả kỳ vọng cuối cùng: khi dự án xong, ta có gì trong tay?

Hình dung sản phẩm cuối như **một bộ đồ nghề, ba bài báo, một bảng kiểm và một luận án**:

1. **Bộ đồ nghề phần mềm `imaihda` phiên bản 1.0.** Cài được từ kho chính thức của R (CRAN), có hướng dẫn từng bước, bản R và Python làm được như nhau. Ai có dữ liệu cũng dùng được, không cần hiểu hết toán bên trong.
2. **Bài báo A (phương pháp):** "So sánh VPC/PCV thô có thể sai; đây là cách kiểm tra và sửa." Dựa trên mô phỏng, viết theo chuẩn báo cáo mô phỏng quốc tế.
3. **Bài báo B (lý thuyết):** "Khi nào có thể tính đúng VPC, khi nào chỉ có thể đưa ra một khoảng." Dành cho tạp chí thống kê.
4. **Bài báo C (thực nghiệm, cốt lõi của luận án):** so sánh bất bình đẳng giao thoa về sức khoẻ ở người trên 50 tuổi giữa nước thu nhập cao và thu nhập trung bình, có dùng bộ đồ nghề để kiểm tra "cân lệch".
5. **Bảng kiểm báo cáo:** danh sách ngắn các việc tác giả nên làm và nên báo cáo khi dùng I-MAIHDA mà nghi ngờ bệnh được phát hiện không đều.
6. **Luận án:** các chương ghép từ các sản phẩm trên.

### Ý tưởng cầu nối quan trọng: dùng xét nghiệm để "đo độ lệch của cân"

Hiện nay mức độ "bệnh bị giấu" (ký hiệu δ) là con số ta phải *giả định*. Nhưng nhiều cohort lão hoá vừa hỏi người tham gia "bác sĩ đã bao giờ nói ông/bà bị tiểu đường chưa?", vừa **xét nghiệm máu** (ví dụ HbA1c) hoặc **đo huyết áp**. So sánh hai nguồn này cho từng nhóm xã hội sẽ cho biết *thật sự* bao nhiêu người bệnh không được chẩn đoán. Như vậy ta **đo được độ lệch của cân** thay vì phải đoán.

Đây là đề xuất; cần kiểm tra cụ thể từng cohort có loại dữ liệu này ở đợt khảo sát nào.

---

## 6. Các bước tiếp theo

Lộ trình gồm 6 giai đoạn. Các giai đoạn có thể chồng lên nhau. Mỗi giai đoạn ghi: làm gì, vì sao, và **khi nào coi là xong**.

### Giai đoạn 0 — Dọn nhà (khoảng 1 tháng)

- **Làm gì:** sửa các lỗi nhỏ đã tìm thấy: đường dẫn cài đặt sai tên repo, tài liệu nhắc tới file không tồn tại, số liệu kiểm thử ghi không thống nhất, bản Python còn dùng công thức cũ, trang hướng dẫn của gói R chưa được đưa vào repo. Viết một lệnh duy nhất chạy toàn bộ kiểm thử **trên máy mình**.
- **Vì sao:** nền móng phải chắc trước khi xây tiếp; người đọc bài báo sẽ cài thử gói.
- **Xong khi:** cài đặt từ GitHub chạy ngay; toàn bộ kiểm thử R và Python đều đạt; README khớp với mã.
- **Việc song song cần làm ngay:** **nộp đơn xin quyền truy cập dữ liệu** các cohort. Thủ tục này thường mất nhiều tháng, nên phải bắt đầu sớm.

### Giai đoạn 1 — Làm cho bản R đầy đủ và phát hành bản 1.0 (khoảng 3 tháng)

- **Làm gì:** chuyển công cụ sửa đồng thời và "khoảng an toàn" từ Python sang R; viết hướng dẫn sử dụng có ví dụ ("thang giả định": khi nào dùng công cụ nào); đưa gói lên CRAN.
- **Vì sao:** người dùng I-MAIHDA chủ yếu dùng R; một gói trên CRAN được tin cậy và trích dẫn nhiều hơn.
- **Xong khi:** gói qua được kiểm tra chính thức của CRAN không lỗi; R và Python cho cùng kết luận trên cùng dữ liệu.

### Giai đoạn 2 — Viết Bài báo A (khoảng 4 tháng)

- **Làm gì:** viết *trước* kế hoạch mô phỏng (mục tiêu, cách sinh dữ liệu, con số cần ước lượng, phương pháp so sánh, cách đo hiệu quả); tăng số lần lặp để kết quả đủ chính xác; so sánh với phương pháp chuẩn đầy đủ; viết bản thảo.
- **Vì sao:** đây là bài dễ hoàn thành nhất vì phần lớn kết quả đã có.
- **Xong khi:** bản thảo được gửi tạp chí.

### Giai đoạn 3 — Viết Bài báo B (khoảng 6 tháng, làm song song với giai đoạn 2)

- **Làm gì:** hoàn tất các chứng minh toán học còn thiếu; làm "khoảng an toàn" cho PCV; cho công cụ chạy với số nhóm bất kỳ (không chỉ 36); rà soát lại tài liệu tham khảo.
- **Vì sao:** đây là đóng góp mới nhất về mặt khoa học. Theo rà soát sơ bộ của dự án, chưa ai xử lý vấn đề "bệnh bị giấu" trong I-MAIHDA; nhận định này sẽ được kiểm chứng lại bằng tổng quan tài liệu có hệ thống.
- **Xong khi:** bản thảo được gửi tạp chí thống kê.

### Giai đoạn 4 — Dùng trên dữ liệu thật (khoảng 6–9 tháng)

- **Làm gì:**
  1. Kiểm tra dữ liệu các cohort có thể ghép được với nhau không (cùng định nghĩa giới, học vấn, tài sản, nơi ở, bệnh).
  2. Đăng ký trước kế hoạch phân tích (để tránh "chọn kết quả đẹp").
  3. Chạy phân tích I-MAIHDA chuẩn, rồi chạy tiếp các bước kiểm tra "cân lệch" theo thang giả định.
  4. Dùng dữ liệu xét nghiệm để đo độ lệch của cân (mục 5).
- **Vì sao:** đây chính là câu hỏi của luận án.
- **Xong khi:** có kết luận so sánh nước giàu – nước thu nhập trung bình, kèm câu trả lời "kết luận này có còn đúng không nếu cân bị lệch đến mức X?".
- **Ứng viên dữ liệu** (cần đối chiếu với exposé): HRS (Mỹ, đã được nêu trong báo cáo QC), SHARE (châu Âu), ELSA (Anh) cho nhóm thu nhập cao; CHARLS (Trung Quốc), LASI (Ấn Độ), MHAS (Mexico), ELSI (Brazil) cho nhóm thu nhập trung bình.

### Giai đoạn 5 — Tổng hợp (khoảng 3 tháng)

- **Làm gì:** hoàn thiện bảng kiểm báo cáo, viết hướng dẫn thực hành, ghép các chương luận án.
- **Xong khi:** luận án hoàn chỉnh; gói phần mềm và bảng kiểm công khai cho cộng đồng.

---

## 7. Lịch trình tổng (dự kiến, bắt đầu từ tháng 10/2026)

| Giai đoạn | Q4/2026 | Q1/2027 | Q2/2027 | Q3/2027 | Q4/2027 | Q1/2028 |
|---|:-:|:-:|:-:|:-:|:-:|:-:|
| 0 · Dọn nhà + nộp đơn xin dữ liệu | ■ | | | | | |
| 1 · Bản R đầy đủ, phát hành 1.0 | ■ | ■ | | | | |
| 2 · Bài báo A | | ■ | ■ | | | |
| 3 · Bài báo B | | ■ | ■ | ■ | | |
| 4 · Dữ liệu thật (Bài báo C) | | | ■ | ■ | ■ | |
| 5 · Tổng hợp, luận án | | | | | ■ | ■ |

Thời gian chỉ mang tính định hướng; phụ thuộc lớn nhất là **thời gian được cấp dữ liệu**.

---

## 8. Những rủi ro và cách phòng

| Rủi ro | Cách phòng |
|---|---|
| Xin dữ liệu lâu | Nộp đơn ngay từ giai đoạn 0; trong lúc chờ, làm bài báo A và B |
| Dữ liệu các nước không ghép được với nhau | Kiểm tra khả năng ghép sớm; dùng bộ dữ liệu đã được hài hoà sẵn nếu có; chấp nhận thu hẹp số nước |
| Không có dữ liệu xét nghiệm để đo độ lệch | Quay về cách "giả định một khoảng hợp lý" (đã có công cụ) và báo cáo "khoảng an toàn" |
| Công cụ cho kết quả sai khi giả định sai | Đã biết rõ vùng nguy hiểm (mục 3); luôn báo cáo kèm phân tích độ nhạy và "khoảng an toàn" |
| Nhận định "chưa ai làm" hoá ra sai | Làm tổng quan tài liệu có hệ thống trước khi viết phần giới thiệu |

---

## 9. Nguyên tắc "không dùng cloud"

- Mọi phân tích **dữ liệu thật** chỉ chạy trên máy cục bộ của người nghiên cứu. Không tải dữ liệu thật lên dịch vụ đám mây, công cụ AI trực tuyến hay máy chủ bên ngoài. Điều này cũng phù hợp với thoả thuận sử dụng dữ liệu của các cohort.
- Kho mã (GitHub) chỉ chứa **mã nguồn, dữ liệu mô phỏng và tài liệu**, không chứa dữ liệu thật.
- Toàn bộ kiểm thử phải chạy được bằng một lệnh trên máy cá nhân, không phụ thuộc dịch vụ bên ngoài.

---

## 10. Nếu chỉ nhớ ba điều

1. **So sánh VPC/PCV thô giữa các nước có thể sai.** Dự án đã chứng minh điều này và đã có công cụ để kiểm tra, sửa, hoặc báo cáo một khoảng trung thực.
2. **Phần còn lại là hoàn thiện và công bố:** dọn nền, làm bản R đầy đủ, viết hai bài báo phương pháp.
3. **Đích cuối là dữ liệu thật:** so sánh nước thu nhập cao với nước thu nhập trung bình ở người trên 50 tuổi, trong đó dùng xét nghiệm máu để đo trực tiếp mức "bệnh bị giấu".

---

*Tài liệu định hướng. Không chứa dữ liệu thật và không đưa ra kết luận thực nghiệm nào về bất kỳ quần thể nào.*
