# Figure 终审报告

**日期**: 2026-07-10  
**范围**: Figure 1–6, Supplementary Figures S1–S4  
**检查项**: 分辨率 / 文件完整性 / 命名一致性 / Caption时态 / 字体 / 颜色 / 图例

---

## 一、文件清单与分辨率

### Main Figures

| 编号 | 文件名 | 像素 | DPI | 状态 |
|------|--------|------|-----|------|
| Figure 1 | Figure1_Study_Design.png | — | 300 | ✅ |
| Figure 2 | Figure2_WGCNA_GSE141198.png (三拼版) | — | 300 | ✅ |
| Figure 2A | Figure2A_dendrogram.png (中间件) | — | 300 | ✅ |
| Figure 2B | Figure2B_GO_Bubble.png (中间件) | — | 300 | ✅ |
| Figure 3 | Figure_ssGSEA_cross_disease.png | — | 300 | ✅ |
| Figure 4 | Figure4_TGS_Validation.png | — | 300 | ⚠️ 内容与caption不匹配 |
| Figure 5A | Figure_TF_TGS_correlation.png | — | 300 | ✅ |
| Figure 5B | Figure_Pathway_TGS_correlation.png | — | 300 | ✅ |
| Figure 6 | Figure6_Mechanistic_Model.png | — | 300 | ✅ |

### Supplementary Figures

| 编号 | 文件名 | 像素 | DPI | 状态 |
|------|--------|------|-----|------|
| Fig S1 | Figure_S2_Direction_Consistency.png | 2700×2400 | 300 | ✅ (命名见问题3) |
| Fig S2A | GSE141198_TGS_KM.png | — | 300 | ✅ |
| Fig S2B | GSE141198_RPL39_KM.png | — | 300 | ✅ |
| Fig S3 | Figure_TF_Fisher_enrichment.png | — | 300 | ✅ |
| **Fig S4A** | **Figure_S1_SoftThreshold_GSE141198.png** | 1500×750 | **150** | 🚨 |
| Fig S4B | Figure_S1_SoftThreshold_GSE57338.png | — | 300 | ✅ |

### 冗余/孤儿文件

| 文件名 | 问题 | 建议 |
|--------|------|------|
| Figure_CrossDisease_DirectionTest.png | 与 Figure_S2_Direction_Consistency.png MD5 完全相同 (8583382be...) | **删除** |
| Figure_CrossDisease_Forest.png | 已被嵌入 Figure4 (panel B)，独立文件未在manuscript引用 | 保留或移入 archive/ |
| Figure2A_dendrogram.png | 已嵌入 Figure2，独立保留 | 保留或移入 archive/ |
| Figure2B_GO_Bubble.png | 已嵌入 Figure2，独立保留 | 保留或移入 archive/ |

---

## 二、CRITICAL 问题

### 🚨 C1: Figure_S1_SoftThreshold_GSE141198.png 分辨率仅 150 DPI

- **现状**: 1500×750 px @ 150 DPI (来源: `WGCNA_GSE141198/soft_threshold.png` 也是 150 DPI)
- **脚本**: `08_generate_remaining_figures.R` L144-147 直接 `file.copy()` 从 WGCNA 目录
- **修复**: 在 `04b_GSE141198_WGCNA_final.R` 中将 `soft_threshold.png` 输出改为 300 DPI，或重建此图为 300 DPI
- **BBA 要求**: ≥300 DPI

### 🚨 C2: Figure 4 内容与 Caption 不匹配

**Caption 声称** (L97-98):
> (A) GSE14520 (n = 221); (B) GSE76427 (n = 115); (C) GSE141198 (n = 148). (D) Forest plot summarizing univariate Cox regression results across cohorts.

**实际文件** (`08_generate_remaining_figures.R` L27-38):
- Panel A: `GSE141198_TGS_KM.png` (仅 GSE141198 一条 KM 曲线)
- Panel B: `Figure_CrossDisease_Forest.png` (森林图)

**缺失**: GSE14520 和 GSE76427 的 Kaplan-Meier 曲线根本没有在这个 figure 里！

**文本交叉引用** (L86):
> "Consistent with the negative validation results observed in GSE14520 (n = 221) and GSE76427 (n = 115) (Figure 4A–B)"

→ 正文明确引用 Figure 4A–B 为 GSE14520 和 GSE76427 的 KM 曲线，但图中不存在。

**选项**:
- A) 补充 GSE14520_TGS_KM.png、GSE76427_TGS_KM.png 并重新组装 Figure 4（4 面板）
- B) 更新 Caption 和正文，改为仅展示 GSE141198 + meta 森林图

---

## 三、MODERATE 问题

### ⚠️ M1: 文件命名与稿件编号不一致

| 实际文件名 | 稿件中的编号 | 混乱点 |
|------------|-------------|--------|
| Figure_S1_SoftThreshold_*.png | **Supplementary Figure S4** | 文件名写 S1，稿件编号 S4 |
| Figure_S2_Direction_Consistency.png | **Supplementary Figure S1** | 文件名写 S2，稿件编号 S1 |

**原因**: `08_generate_remaining_figures.R` 中的注释按生成顺序编号（先做的软阈值→S1，后做的方向一致性→S2），但稿件排版时重新排序（按出现顺序：方向一致性在 Results 2.2 先出现→S1，软阈值在 Methods 4.2 后出现→S4）。

**建议**: 重命名文件以匹配稿件编号：
- `Figure_S1_SoftThreshold_GSE141198.png` → `Figure_S4A_SoftThreshold_GSE141198.png`
- `Figure_S1_SoftThreshold_GSE57338.png` → `Figure_S4B_SoftThreshold_GSE57338.png`
- `Figure_S2_Direction_Consistency.png` → `Figure_S1_Direction_Consistency.png`

### ⚠️ M2: 字体大小不统一

| 图 | base_size | 备注 |
|----|-----------|------|
| Figure 2B (GO bubble) | 11 | 最小 |
| Figure 4 forest (Figure_CrossDisease_Forest) | 12 | |
| Figure 5A (TF-TGS) | 12 | |
| Figure 5B (Pathway-TGS) | 12 | |
| Figure S3 (Fisher) | 12 | |
| Figure 3 (ssGSEA scatter) | 13 | |
| Figure S2A/B (KM curves) | 13 | |
| Figure_CrossDisease_Forest | 13 | |
| Figure S1 (Direction Consistency) | 14 | 最大 |

**范围**: 11–14。BBA 无硬性规定，但建议统一为 12 或 13。

### ⚠️ M3: 未显式指定字体族

所有 ggplot2 图使用 `theme_minimal()` 但未设置 `base_family`。Windows 下默认渲染为 Arial（sans-serif）。对于 BBA 投稿通常可接受，但建议显式设置：
```r
theme_minimal(base_family = "sans")  # = Arial on Windows
```
如果要在 PDF 中嵌入字体，需 `extrafont` 包加载 Arial。

### ⚠️ M4: Figure 2B GO 气泡颜色映射

使用 `scale_color_gradient(low = "red", high = "steelblue")`，p.adjust 从红渐变到蓝 → 中间值呈紫色，与稿件中其他图的色系不一致。建议改为：
```r
scale_color_gradient(low = "#E41A1C", high = "#BDBDBD")  # 红→灰
```

### ⚠️ M5: Figure 3 Caption 象限描述不全

> "Quadrant I (upper right, HCC↑ HF↑): 2 pathways; Quadrant IV (lower right, HCC↑ HF↓): 31 of 33 translation pathways."

只描述了 I 和 IV 象限，未提及 II（左上, HCC↓ HF↑）和 III（左下, HCC↓ HF↓）。如果这两个象限数据不重要或为空，建议加一句如 "No pathways fell in Quadrants II or III"。

---

## 四、MINOR 问题

### ℹ️ I1: Figure 2A 使用 Base R 图形

- `plotDendroAndColors()` 是 WGCNA 自带 base R 函数，输出 PNG 设备
- 所有其他图使用 ggplot2 → 渲染风格不完全一致（base R 标题字体略大、间距不同）
- **可接受**（WGCNA dendrogram 无 ggplot 替代），无需修改

### ℹ️ I2: Figure 6 颜色编码

- HCC 侧用红色 (`#B71C1C`)，HF 侧用蓝色 (`#0D47A1`)
- 跨疾病分析 Figure 3 同样用红色标记翻译通路、灰色标记非翻译通路
- 色系一致 ✅

### ℹ️ I3: Caption 时态

全部统一使用现在时（descriptive present tense），格式为 `**Figure X. Title.** Description.`。✅ 一致。

---

## 五、Figure-Caption 对应总表

| 稿件编号 | 文件 | Caption位置 | 正文引用 | 状态 |
|----------|------|------------|---------|------|
| Fig 1 | Figure1_Study_Design.png | L38 | — | ✅ |
| Fig 2 | Figure2_WGCNA_GSE141198.png | L59 | L48 | ✅ (A=树状图, B=GO, C=热图) |
| Fig 3 | Figure_ssGSEA_cross_disease.png | L80 | L65-68 | ⚠️ 象限描述不全 (M5) |
| Fig 4 | Figure4_TGS_Validation.png | L97 | L86 | 🚨 缺GSE14520/GSE76427 KM (C2) |
| Fig 5 | Figure_TF_TGS_correlation.png + Figure_Pathway_TGS_correlation.png | L158 | L123-125 | ✅ |
| Fig 6 | Figure6_Mechanistic_Model.png | L198 | — | ✅ |
| Fig S1 | Figure_S2_Direction_Consistency.png | L169 | L69 | ⚠️ 文件名S2→应为S1 (M1) |
| Fig S2 | GSE141198_TGS_KM.png + GSE141198_RPL39_KM.png | L117 | L86 | ✅ |
| Fig S3 | Figure_TF_Fisher_enrichment.png | L138 | L127 | ✅ |
| Fig S4 | Figure_S1_SoftThreshold_GSE*.png | L230 | — | 🚨 GSE141198 150 DPI (C1) + 文件名S1→应为S4 (M1) |

---

## 六、修复优先级

| 优先级 | 问题 | 修复方式 |
|--------|------|---------|
| 🔴 P0 | C1: S4A 分辨率 150 DPI | 重新生成 GSE141198 软阈值图为 300 DPI |
| 🔴 P0 | C2: Figure 4 内容不完整 | 补充 GSE14520/GSE76427 KM 曲线 OR 修正 caption |
| 🟡 P1 | M1: 文件命名 | 统一重命名为与稿件编号一致 |
| 🟡 P1 | M5: Fig 3 Caption 象限 | 补充描述 or 注明无数据 |
| 🟢 P2 | M2: 字体大小 | 统一 base_size |
| 🟢 P2 | M3: 字体族 | 显式设置 base_family |
| 🟢 P2 | M4: GO 气泡颜色 | 改用单一色渐变 |
| ⚪ P3 | 冗余文件清理 | 删除 Duplicate Figure_CrossDisease_DirectionTest.png |
