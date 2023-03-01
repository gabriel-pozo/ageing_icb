#Load 10x data

MD01024 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352886_MD01-024_tumor_1/MD01-024_tumor_1/")
MD01024 <- CreateSeuratObject(counts = MD01024, project = 'MD01024', min.cells = 3, min.features = 200)

MD01010 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352888_MD01-010_tumor_1/MD01-010_tumor_1/")
MD01010 <- CreateSeuratObject(counts = MD01010, project = 'MD01010', min.cells = 3, min.features = 200)

MD01004 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352889_MD01-004_tumor_1/MD01-004_tumor_1/")
MD01004 <- CreateSeuratObject(counts = MD01004, project = 'MD01004', min.cells = 3, min.features = 200)

MD01004 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352889_MD01-004_tumor_1/MD01-004_tumor_1/")
MD01004 <- CreateSeuratObject(counts = MD01004, project = 'MD01004', min.cells = 3, min.features = 200)

MD043011M1 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352893_MD043-011_mettumor_1/MD043-011_mettumor_1/")
MD043011M1 <- CreateSeuratObject(counts = MD043011M1, project = 'MD043011M1', min.cells = 3, min.features = 200)

MD043011M2 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352895_MD043-011_mettumor_2/MD043-011_mettumor_2/")
MD043011M2 <- CreateSeuratObject(counts = MD043011M2, project = 'MD043011M2', min.cells = 3, min.features = 200)

MD043011M3 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352896_MD043-011_mettumor_3/MD043-011_mettumor_3/")
MD043011M3 <- CreateSeuratObject(counts = MD043011M3, project = 'MD043011M3', min.cells = 3, min.features = 200)

MD043011M4 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352897_MD043-011_mettumor_4/MD043-011_mettumor_4/")
MD043011M4 <- CreateSeuratObject(counts = MD043011M4, project = 'MD043011M4', min.cells = 3, min.features = 200)

MD043011M5 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352898_MD043-011_mettumor_5/MD043-011_mettumor_5/")
MD043011M5 <- CreateSeuratObject(counts = MD043011M5, project = 'MD043011M5', min.cells = 3, min.features = 200)

MD043011T2 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352899_MD043-011_tumor_2/MD043-011_tumor_2/")
MD043011T2 <- CreateSeuratObject(counts = MD043011T2, project = 'MD043011T2', min.cells = 3, min.features = 200)

MD043011T3 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352900_MD043-011_tumor_3/MD043-011_tumor_3/")
MD043011T3 <- CreateSeuratObject(counts = MD043011T3, project = 'MD043011T3', min.cells = 3, min.features = 200)

MD043011T4 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352901_MD043-011_tumor_4/MD043-011_tumor_4/")
MD043011T4 <- CreateSeuratObject(counts = MD043011T4, project = 'MD043011T4', min.cells = 3, min.features = 200)

MD043011T5 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352902_MD043-011_tumor_5/MD043-011_tumor_5/")
MD043011T5 <- CreateSeuratObject(counts = MD043011T5, project = 'MD043011T5', min.cells = 3, min.features = 200)

MD01019T1 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352910_MD01-019_tumor_1/MD01-019_tumor_1/")
MD01019T1 <- CreateSeuratObject(counts = MD01019T1, project = 'MD01019T1', min.cells = 3, min.features = 200)

MD01019T10 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352912_MD01-019_tumor_10/MD01-019_tumor_10/")
MD01019T10 <- CreateSeuratObject(counts = MD01019T10, project = 'MD01019T10', min.cells = 3, min.features = 200)

MD01019T11 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352913_MD01-019_tumor_11/MD01-019_tumor_11/")
MD01019T11 <- CreateSeuratObject(counts = MD01019T11, project = 'MD01019T11', min.cells = 3, min.features = 200)

MD01019T13 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352914_MD01-019_tumor_13/MD01-019_tumor_13/")
MD01019T13 <- CreateSeuratObject(counts = MD01019T13, project = 'MD01019T13', min.cells = 3, min.features = 200)

MD01019T15 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352915_MD01-019_tumor_15/MD01-019_tumor_15/")
MD01019T15 <- CreateSeuratObject(counts = MD01019T15, project = 'MD01019T15', min.cells = 3, min.features = 200)

MD01019T17 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352916_MD01-019_tumor_17/MD01-019_tumor_17/")
MD01019T17 <- CreateSeuratObject(counts = MD01019T17, project = 'MD01019T17', min.cells = 3, min.features = 200)

MD01019T3 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352917_MD01-019_tumor_3/MD01-019_tumor_3/")
MD01019T3 <- CreateSeuratObject(counts = MD01019T3, project = 'MD01019T3', min.cells = 3, min.features = 200)

MD01019T6 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352918_MD01-019_tumor_6/MD01-019_tumor_6/")
MD01019T6 <- CreateSeuratObject(counts = MD01019T6, project = 'MD01019T6', min.cells = 3, min.features = 200)

MD01019T7 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352919_MD01-019_tumor_7/MD01-019_tumor_7/")
MD01019T7 <- CreateSeuratObject(counts = MD01019T7, project = 'MD01019T7', min.cells = 3, min.features = 200)

MD01019T8 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352920_MD01-019_tumor_8/MD01-019_tumor_8/")
MD01019T8 <- CreateSeuratObject(counts = MD01019T8, project = 'MD01019T8', min.cells = 3, min.features = 200)

MD01019T25 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352922_MD01-019_tumor_25/MD01-019_tumor_25/")
MD01019T25 <- CreateSeuratObject(counts = MD01019T25, project = 'MD01019T25', min.cells = 3, min.features = 200)

NY016007T1 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352924_NY016-007_tumor_1/NY016-007_tumor_1/")
NY016007T1 <- CreateSeuratObject(counts = NY016007T1, project = 'NY016007T1', min.cells = 3, min.features = 200)

NY016007T2 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352925_NY016-007_tumor_2/NY016-007_tumor_2/")
NY016007T2 <- CreateSeuratObject(counts = NY016007T2, project = 'NY016007T2', min.cells = 3, min.features = 200)

NY016007T3 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352926_NY016-007_tumor_3/NY016-007_tumor_3/")
NY016007T3 <- CreateSeuratObject(counts = NY016007T3, project = 'NY016007T3', min.cells = 3, min.features = 200)

NY016014T1 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352927_NY016-014_tumor_1/NY016-014_tumor_1/")
NY016014T1 <- CreateSeuratObject(counts = NY016014T1, project = 'NY016014T1', min.cells = 3, min.features = 200)

NY016014T2 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352928_NY016-014_tumor_2/NY016-014_tumor_2/")
NY016014T2 <- CreateSeuratObject(counts = NY016014T2, project = 'NY016014T2', min.cells = 3, min.features = 200)

NY016014T3 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352929_NY016-014_tumor_3/NY016-014_tumor_3/")
NY016014T3 <- CreateSeuratObject(counts = NY016014T3, project = 'NY016014T3', min.cells = 3, min.features = 200)

NY016015T1 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352931_NY016-015_tumor_1/NY016-015_tumor_1/")
NY016015T1 <- CreateSeuratObject(counts = NY016015T1, project = 'NY016015T1', min.cells = 3, min.features = 200)

NY016015T2 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352932_NY016-015_tumor_2/NY016-015_tumor_2/")
NY016015T2 <- CreateSeuratObject(counts = NY016015T2, project = 'NY016015T2', min.cells = 3, min.features = 200)

NY016015T3 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352933_NY016-015_tumor_3/NY016-015_tumor_3/")
NY016015T3 <- CreateSeuratObject(counts = NY016015T3, project = 'NY016015T3', min.cells = 3, min.features = 200)

NY016015T4 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352934_NY016-015_tumor_4/NY016-015_tumor_4/")
NY016015T4 <- CreateSeuratObject(counts = NY016015T4, project = 'NY016015T4', min.cells = 3, min.features = 200)

NY016021 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352936_NY016-021_tumor_1/NY016-021_tumor_1/")
NY016021 <- CreateSeuratObject(counts = NY016021, project = 'NY016021', min.cells = 3, min.features = 200)

NY016016 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352935_NY016-016_normal_1/NY016-016_normal_1/")
NY016016 <- CreateSeuratObject(counts = NY016016, project = 'NY016016', min.cells = 3, min.features = 200)

NY016022T1 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352938_NY016-022_tumor_1/NY016-022_tumor_1/")
NY016022T1 <- CreateSeuratObject(counts = NY016022T1, project = 'NY016022T1', min.cells = 3, min.features = 200)

NY016022T2 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352939_NY016-022_tumor_2/NY016-022_tumor_2/")
NY016022T2 <- CreateSeuratObject(counts = NY016022T2, project = 'NY016022T2', min.cells = 3, min.features = 200)

NY016022T3 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352940_NY016-022_tumor_3/NY016-022_tumor_3/")
NY016022T3 <- CreateSeuratObject(counts = NY016022T3, project = 'NY016022T3', min.cells = 3, min.features = 200)

NY016022T4 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352941_NY016-022_tumor_4/NY016-022_tumor_4/")
NY016022T4 <- CreateSeuratObject(counts = NY016022T4, project = 'NY016022T4', min.cells = 3, min.features = 200)

NY016025T1 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352943_NY016-025_tumor_1/NY016-025_tumor_1/")
NY016025T1 <- CreateSeuratObject(counts = NY016025T1, project = 'NY016025T1', min.cells = 3, min.features = 200)

NY016025T2 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352944_NY016-025_tumor_2/NY016-025_tumor_2/")
NY016025T2 <- CreateSeuratObject(counts = NY016025T2, project = 'NY016025T2', min.cells = 3, min.features = 200)

NY016025T3 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352945_NY016-025_tumor_3/NY016-025_tumor_3/")
NY016025T3 <- CreateSeuratObject(counts = NY016025T3, project = 'NY016025T3', min.cells = 3, min.features = 200)

NY016025T4 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352946_NY016-025_tumor_4/NY016-025_tumor_4/")
NY016025T4 <- CreateSeuratObject(counts = NY016025T4, project = 'NY016025T4', min.cells = 3, min.features = 200)

MD043008 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352947_MD043-008_tumor_1/MD043-008_tumor_1/")
MD043008 <- CreateSeuratObject(counts = MD043008, project = 'MD043008', min.cells = 3, min.features = 200)

MD043003T1 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352949_MD043-003_tumor_1/MD043-003_tumor_1/")
MD043003T1 <- CreateSeuratObject(counts = MD043003T1, project = 'MD043003T1', min.cells = 3, min.features = 200)

MD043003T2 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352951_MD043-003_tumor_2/MD043-003_tumor_2/")
MD043003T2 <- CreateSeuratObject(counts = MD043003T2, project = 'MD043003T2', min.cells = 3, min.features = 200)

MD043003T3 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352952_MD043-003_tumor_3/MD043-003_tumor_3/")
MD043003T3 <- CreateSeuratObject(counts = MD043003T3, project = 'MD043003T3', min.cells = 3, min.features = 200)

MD043003T4 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352953_MD043-003_tumor_4/MD043-003_tumor_4/")
MD043003T4 <- CreateSeuratObject(counts = MD043003T4, project = 'MD043003T4', min.cells = 3, min.features = 200)

MD01005T2 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352960_MD01-005_tumor_2/MD01-005_tumor_2/")
MD01005T2 <- CreateSeuratObject(counts = MD01005T2, project = 'MD01005T2', min.cells = 3, min.features = 200)

MD01005T3 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352961_MD01-005_tumor_3/MD01-005_tumor_3/")
MD01005T3 <- CreateSeuratObject(counts = MD01005T3, project = 'MD01005T3', min.cells = 3, min.features = 200)

MD01005T4 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352962_MD01-005_tumor_4/MD01-005_tumor_4/")
MD01005T4 <- CreateSeuratObject(counts = MD01005T4, project = 'MD01005T4', min.cells = 3, min.features = 200)

MD01005T5 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352963_MD01-005_tumor_5/MD01-005_tumor_5/")
MD01005T5 <- CreateSeuratObject(counts = MD01005T5, project = 'MD01005T5', min.cells = 3, min.features = 200)

MD01005T6 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352964_MD01-005_tumor_6/MD01-005_tumor_6/")
MD01005T6 <- CreateSeuratObject(counts = MD01005T6, project = 'MD01005T6', min.cells = 3, min.features = 200)

MD01005T7 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352965_MD01-005_tumor_7/MD01-005_tumor_7/")
MD01005T7 <- CreateSeuratObject(counts = MD01005T7, project = 'MD01005T7', min.cells = 3, min.features = 200)

MD01005T8 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352966_MD01-005_tumor_8/MD01-005_tumor_8/")
MD01005T8 <- CreateSeuratObject(counts = MD01005T8, project = 'MD01005T8', min.cells = 3, min.features = 200)

MD01005T9 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352967_MD01-005_tumor_9/MD01-005_tumor_9/")
MD01005T9 <- CreateSeuratObject(counts = MD01005T9, project = 'MD01005T9', min.cells = 3, min.features = 200)

MD043006 <- Read10X(data.dir = "D:/Scanpy/Lung/GSE173351_RAW/GSM5352973_MD043-006_tumor_1/MD043-006_tumor_1/")
MD043006 <- CreateSeuratObject(counts = MD043006, project = 'MD043006', min.cells = 3, min.features = 200)

#Add metadata information

MD01004$ICB <- "aPD1"
MD01004$Sex <- "M"
MD01004$Age <- "67"
MD01004$ICB_response <- "non-MPR"
MD01004$Stage <- "IV"
MD01004$sample_id <- "MD01004"
MD01004$patient_id <- "MD01004"

MD01010$ICB <- "aPD1"
MD01010$Sex <- "F"
MD01010$Age <- "78"
MD01010$ICB_response <- "MPR"
MD01010$Stage <- "III"
MD01010$sample_id <- "MD01010"
MD01010$patient_id <- "MD01010"

MD01024$ICB <- "aPD1"
MD01024$Sex <- "F"
MD01024$Age <- "70"
MD01024$ICB_response <- "non-MPR"
MD01024$Stage <- "I"
MD01024$sample_id <- "MD01024"
MD01024$patient_id <- "MD01024"

MD043006$ICB <- "aPD1"
MD043006$Sex <- "M"
MD043006$Age <- "69"
MD043006$ICB_response <- "non-MPR"
MD043006$Stage <- "II"
MD043006$sample_id <- "MD043006"
MD043006$patient_id <- "MD043006"

MD043008$ICB <- "aPD1"
MD043008$Sex <- "F"
MD043008$Age <- "72"
MD043008$ICB_response <- "MPR"
MD043008$Stage <- "I"
MD043008$sample_id <- "MD043008"
MD043008$patient_id <- "MD043008"

MD01005T2$ICB <- "aPD1"
MD01005T2$Sex <- "M"
MD01005T2$Age <- "61"
MD01005T2$ICB_response <- "MPR"
MD01005T2$Stage <- "III"
MD01005T2$sample_id <- "MD01005T2"
MD01005T2$patient_id <- "MD01005"

MD01005T3$ICB <- "aPD1"
MD01005T3$Sex <- "M"
MD01005T3$Age <- "61"
MD01005T3$ICB_response <- "MPR"
MD01005T3$Stage <- "III"
MD01005T3$sample_id <- "MD01005T3"
MD01005T3$patient_id <- "MD01005"

MD01005T4$ICB <- "aPD1"
MD01005T4$Sex <- "M"
MD01005T4$Age <- "61"
MD01005T4$ICB_response <- "MPR"
MD01005T4$Stage <- "III"
MD01005T4$sample_id <- "MD01005T4"
MD01005T4$patient_id <- "MD01005"

MD01005T5$ICB <- "aPD1"
MD01005T5$Sex <- "M"
MD01005T5$Age <- "61"
MD01005T5$ICB_response <- "MPR"
MD01005T5$Stage <- "III"
MD01005T5$sample_id <- "MD01005T5"
MD01005T5$patient_id <- "MD01005"

MD01005T6$ICB <- "aPD1"
MD01005T6$Sex <- "M"
MD01005T6$Age <- "61"
MD01005T6$ICB_response <- "MPR"
MD01005T6$Stage <- "III"
MD01005T6$sample_id <- "MD01005T6"
MD01005T6$patient_id <- "MD01005"

MD01005T7$ICB <- "aPD1"
MD01005T7$Sex <- "M"
MD01005T7$Age <- "61"
MD01005T7$ICB_response <- "MPR"
MD01005T7$Stage <- "III"
MD01005T7$sample_id <- "MD01005T7"
MD01005T7$patient_id <- "MD01005"

MD01005T8$ICB <- "aPD1"
MD01005T8$Sex <- "M"
MD01005T8$Age <- "61"
MD01005T8$ICB_response <- "MPR"
MD01005T8$Stage <- "III"
MD01005T8$sample_id <- "MD01005T8"
MD01005T8$patient_id <- "MD01005"

MD01005T9$ICB <- "aPD1"
MD01005T9$Sex <- "M"
MD01005T9$Age <- "61"
MD01005T9$ICB_response <- "MPR"
MD01005T9$Stage <- "III"
MD01005T9$sample_id <- "MD01005T9"
MD01005T9$patient_id <- "MD01005"

MD01019T1$ICB <- "aPD1"
MD01019T1$Sex <- "M"
MD01019T1$Age <- "70"
MD01019T1$ICB_response <- "non-MPR"
MD01019T1$Stage <- "II"
MD01019T1$sample_id <- "MD01019T1"
MD01019T1$patient_id <- "MD01019"

MD01019T3$ICB <- "aPD1"
MD01019T3$Sex <- "M"
MD01019T3$Age <- "70"
MD01019T3$ICB_response <- "non-MPR"
MD01019T3$Stage <- "II"
MD01019T3$sample_id <- "MD01019T3"
MD01019T3$patient_id <- "MD01019"

MD01019T6$ICB <- "aPD1"
MD01019T6$Sex <- "M"
MD01019T6$Age <- "70"
MD01019T6$ICB_response <- "non-MPR"
MD01019T6$Stage <- "II"
MD01019T6$sample_id <- "MD01019T6"
MD01019T6$patient_id <- "MD01019"

MD01019T7$ICB <- "aPD1"
MD01019T7$Sex <- "M"
MD01019T7$Age <- "70"
MD01019T7$ICB_response <- "non-MPR"
MD01019T7$Stage <- "II"
MD01019T7$sample_id <- "MD01019T7"
MD01019T7$patient_id <- "MD01019"

MD01019T8$ICB <- "aPD1"
MD01019T8$Sex <- "M"
MD01019T8$Age <- "70"
MD01019T8$ICB_response <- "non-MPR"
MD01019T8$Stage <- "II"
MD01019T8$sample_id <- "MD01019T8"
MD01019T8$patient_id <- "MD01019"

MD01019T10$ICB <- "aPD1"
MD01019T10$Sex <- "M"
MD01019T10$Age <- "70"
MD01019T10$ICB_response <- "non-MPR"
MD01019T10$Stage <- "II"
MD01019T10$sample_id <- "MD01019T10"
MD01019T10$patient_id <- "MD01019"

MD01019T11$ICB <- "aPD1"
MD01019T11$Sex <- "M"
MD01019T11$Age <- "70"
MD01019T11$ICB_response <- "non-MPR"
MD01019T11$Stage <- "II"
MD01019T11$sample_id <- "MD01019T11"
MD01019T11$patient_id <- "MD01019"

MD01019T13$ICB <- "aPD1"
MD01019T13$Sex <- "M"
MD01019T13$Age <- "70"
MD01019T13$ICB_response <- "non-MPR"
MD01019T13$Stage <- "II"
MD01019T13$sample_id <- "MD01019T13"
MD01019T13$patient_id <- "MD01019"

MD01019T15$ICB <- "aPD1"
MD01019T15$Sex <- "M"
MD01019T15$Age <- "70"
MD01019T15$ICB_response <- "non-MPR"
MD01019T15$Stage <- "II"
MD01019T15$sample_id <- "MD01019T15"
MD01019T15$patient_id <- "MD01019"

MD01019T17$ICB <- "aPD1"
MD01019T17$Sex <- "M"
MD01019T17$Age <- "70"
MD01019T17$ICB_response <- "non-MPR"
MD01019T17$Stage <- "II"
MD01019T17$sample_id <- "MD01019T17"
MD01019T17$patient_id <- "MD01019"

MD01019T25$ICB <- "aPD1"
MD01019T25$Sex <- "M"
MD01019T25$Age <- "70"
MD01019T25$ICB_response <- "non-MPR"
MD01019T25$Stage <- "II"
MD01019T25$sample_id <- "MD01019T25"
MD01019T25$patient_id <- "MD01019"

MD043003T1$ICB <- "aPD1"
MD043003T1$Sex <- "M"
MD043003T1$Age <- "62"
MD043003T1$ICB_response <- "MPR"
MD043003T1$Stage <- "II"
MD043003T1$sample_id <- "MD043003T1"
MD043003T1$patient_id <- "MD043003"

MD043003T2$ICB <- "aPD1"
MD043003T2$Sex <- "M"
MD043003T2$Age <- "62"
MD043003T2$ICB_response <- "MPR"
MD043003T2$Stage <- "II"
MD043003T2$sample_id <- "MD043003T2"
MD043003T2$patient_id <- "MD043003"

MD043003T3$ICB <- "aPD1"
MD043003T3$Sex <- "M"
MD043003T3$Age <- "62"
MD043003T3$ICB_response <- "MPR"
MD043003T3$Stage <- "II"
MD043003T3$sample_id <- "MD043003T3"
MD043003T3$patient_id <- "MD043003"

MD043003T4$ICB <- "aPD1"
MD043003T4$Sex <- "M"
MD043003T4$Age <- "62"
MD043003T4$ICB_response <- "MPR"
MD043003T4$Stage <- "II"
MD043003T4$sample_id <- "MD043003T4"
MD043003T4$patient_id <- "MD043003"

MD043003T5$ICB <- "aPD1"
MD043003T5$Sex <- "M"
MD043003T5$Age <- "62"
MD043003T5$ICB_response <- "MPR"
MD043003T5$Stage <- "II"
MD043003T5$sample_id <- "MD043003T5"
MD043003T5$patient_id <- "MD043003"

MD043011T2$ICB <- "aPD1"
MD043011T2$Sex <- "M"
MD043011T2$Age <- "55"
MD043011T2$ICB_response <- "non-MPR"
MD043011T2$Stage <- "II"
MD043011T2$sample_id <- "MD043011T2"
MD043011T2$patient_id <- "MD043011"

MD043011T3$ICB <- "aPD1"
MD043011T3$Sex <- "M"
MD043011T3$Age <- "55"
MD043011T3$ICB_response <- "non-MPR"
MD043011T3$Stage <- "II"
MD043011T3$sample_id <- "MD043011T3"
MD043011T3$patient_id <- "MD043011"

MD043011T4$ICB <- "aPD1"
MD043011T4$Sex <- "M"
MD043011T4$Age <- "55"
MD043011T4$ICB_response <- "non-MPR"
MD043011T4$Stage <- "II"
MD043011T4$sample_id <- "MD043011T4"
MD043011T4$patient_id <- "MD043011"

MD043011T5$ICB <- "aPD1"
MD043011T5$Sex <- "M"
MD043011T5$Age <- "55"
MD043011T5$ICB_response <- "non-MPR"
MD043011T5$Stage <- "II"
MD043011T5$sample_id <- "MD043011T5"
MD043011T5$patient_id <- "MD043011"

MD043011M1$ICB <- "aPD1"
MD043011M1$Sex <- "M"
MD043011M1$Age <- "55"
MD043011M1$ICB_response <- "non-MPR"
MD043011M1$Stage <- "II"
MD043011M1$sample_id <- "MD043011M1"
MD043011M1$patient_id <- "MD043011"

MD043011M2$ICB <- "aPD1"
MD043011M2$Sex <- "M"
MD043011M2$Age <- "55"
MD043011M2$ICB_response <- "non-MPR"
MD043011M2$Stage <- "II"
MD043011M2$sample_id <- "MD043011M2"
MD043011M2$patient_id <- "MD043011"

MD043011M3$ICB <- "aPD1"
MD043011M3$Sex <- "M"
MD043011M3$Age <- "55"
MD043011M3$ICB_response <- "non-MPR"
MD043011M3$Stage <- "II"
MD043011M3$sample_id <- "MD043011M3"
MD043011M3$patient_id <- "MD043011"

MD043011M4$ICB <- "aPD1"
MD043011M4$Sex <- "M"
MD043011M4$Age <- "55"
MD043011M4$ICB_response <- "non-MPR"
MD043011M4$Stage <- "II"
MD043011M4$sample_id <- "MD043011M4"
MD043011M4$patient_id <- "MD043011"

MD043011M5$ICB <- "aPD1"
MD043011M5$Sex <- "M"
MD043011M5$Age <- "55"
MD043011M5$ICB_response <- "non-MPR"
MD043011M5$Stage <- "II"
MD043011M5$sample_id <- "MD043011M5"
MD043011M5$patient_id <- "MD043011"

#Merge objects and create combined seurat obj

MD <- merge(MD01004, y = c(MD01010, MD01024, MD043006, MD043008, MD01005T2, MD01005T3, MD01005T4,
                           MD01005T5, MD01005T6, MD01005T7, MD01005T8, MD01005T9, MD01019T1, MD01019T3,
                           MD01019T6, MD01019T7, MD01019T8, MD01019T10, MD01019T11, MD01019T13,
                           MD01019T15, MD01019T17, MD01019T25, MD043003T1, MD043003T2, MD043003T3,
                           MD043003T4, MD043003T5, MD043011T2, MD043011T3, MD043011T4, MD043011T5,
                           MD043011M1, MD043011M2, MD043011M3, MD043011M4, MD043011M5), project = "MD",
            add.cell.ids = c("MD01004", "MD01010", "MD01024", "MD043006", "MD043008", "MD01005T2",
                             "MD01005T3", "MD01005T4", "MD01005T5", "MD01005T6", "MD01005T7", "MD01005T8",
                             "MD01005T9", "MD01019T1", "MD01019T3", "MD01019T6", "MD01019T7", "MD01019T8",
                             "MD01019T10", "MD01019T11", "MD01019T13", "MD01019T15", "MD01019T17", "MD01019T25",
                             "MD043003T1", "MD043003T2", "MD043003T3", "MD043003T4", "MD043003T5", "MD043011T2",
                             "MD043011T3", "MD043011T4", "MD043011T5", "MD043011M1", "MD043011M2", "MD043011M3",
                             "MD043011M4", "MD043011M5"))

#QC estimation

MD[["percent.mt"]] <- PercentageFeatureSet(MD, pattern = "^MT-")
MD[["percent.hb"]] <- PercentageFeatureSet(MD, pattern = "^HBA|^HBB")
MD[["percent.rp"]] <- PercentageFeatureSet(MD, pattern = "^RPS|^RPL")


nFeature_lower <- 500
nFeature_upper <- 10000
nCount_lower <- 1000
nCount_upper <- 100000
percent.mt_lower <- 0
percent.mt_upper <- 30
percent.hb_lower <- 0
percent.hb_upper <- 5

#Filter low quality cells
MD <- subset(MD, subset = nFeature_RNA > nFeature_lower & nFeature_RNA < nFeature_upper & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper & percent.mt < percent.mt_upper & percent.hb < percent.hb_upper)

#Prepare metadata and counts matrix for Scanpy
MD$barcode <- colnames(MD)
write.csv(MD@meta.data, file='D:Scanpy/md_anderson/metadata.csv', quote=F, row.names=F)

counts_matrix <- GetAssayData(MD, assay='RNA', slot='counts')
writeMM(counts_matrix, file='D:Scanpy/md_anderson/counts.mtx')

write.table(
  data.frame('gene'=rownames(counts_matrix)),file='D:Scanpy/md_anderson/genes.csv',
  quote=F,row.names=F,col.names=F
)



#Add metadata information

NY016007T1$ICB <- "aPD1"
NY016007T1$Sex <- "F"
NY016007T1$Age <- "68"
NY016007T1$ICB_response <- "non-MPR"
NY016007T1$Stage <- "II"
NY016007T1$sample_id <- "NY016007T1"
NY016007T1$patient_id <- "NY016007"

NY016007T2$ICB <- "aPD1"
NY016007T2$Sex <- "F"
NY016007T2$Age <- "68"
NY016007T2$ICB_response <- "non-MPR"
NY016007T2$Stage <- "II"
NY016007T2$sample_id <- "NY016007T2"
NY016007T2$patient_id <- "NY016007"

NY016007T3$ICB <- "aPD1"
NY016007T3$Sex <- "F"
NY016007T3$Age <- "68"
NY016007T3$ICB_response <- "non-MPR"
NY016007T3$Stage <- "II"
NY016007T3$sample_id <- "NY016007T3"
NY016007T3$patient_id <- "NY016007"

NY016014T1$ICB <- "aPD1"
NY016014T1$Sex <- "F"
NY016014T1$Age <- "58"
NY016014T1$ICB_response <- "non-MPR"
NY016014T1$Stage <- "II"
NY016014T1$sample_id <- "NY016014T1"
NY016014T1$patient_id <- "NY016014"

NY016014T2$ICB <- "aPD1"
NY016014T2$Sex <- "F"
NY016014T2$Age <- "58"
NY016014T2$ICB_response <- "non-MPR"
NY016014T2$Stage <- "II"
NY016014T2$sample_id <- "NY016014T2"
NY016014T2$patient_id <- "NY016014"

NY016014T3$ICB <- "aPD1"
NY016014T3$Sex <- "F"
NY016014T3$Age <- "58"
NY016014T3$ICB_response <- "non-MPR"
NY016014T3$Stage <- "II"
NY016014T3$sample_id <- "NY016014T3"
NY016014T3$patient_id <- "NY016014"

NY016016$ICB <- "aPD1"
NY016016$Sex <- "F"
NY016016$Age <- "79"
NY016016$ICB_response <- "MPR"
NY016016$Stage <- "I"
NY016016$sample_id <- "NY016016"
NY016016$patient_id <- "NY016016"

NY016015T1$ICB <- "aPD1"
NY016015T1$Sex <- "F"
NY016015T1$Age <- "58"
NY016015T1$ICB_response <- "non-MPR"
NY016015T1$Stage <- "II"
NY016015T1$sample_id <- "NY016015T1"
NY016015T1$patient_id <- "NY016015"

NY016015T2$ICB <- "aPD1"
NY016015T2$Sex <- "F"
NY016015T2$Age <- "58"
NY016015T2$ICB_response <- "non-MPR"
NY016015T2$Stage <- "II"
NY016015T2$sample_id <- "NY016015T2"
NY016015T2$patient_id <- "NY016015"

NY016015T3$ICB <- "aPD1"
NY016015T3$Sex <- "F"
NY016015T3$Age <- "58"
NY016015T3$ICB_response <- "non-MPR"
NY016015T3$Stage <- "II"
NY016015T3$sample_id <- "NY016015T3"
NY016015T3$patient_id <- "NY016015"

NY016015T4$ICB <- "aPD1"
NY016015T4$Sex <- "F"
NY016015T4$Age <- "58"
NY016015T4$ICB_response <- "non-MPR"
NY016015T4$Stage <- "II"
NY016015T4$sample_id <- "NY016015T4"
NY016015T4$patient_id <- "NY016015"

NY016021$ICB <- "aPD1"
NY016021$Sex <- "M"
NY016021$Age <- "74"
NY016021$ICB_response <- "non-MPR"
NY016021$Stage <- "III"
NY016021$sample_id <- "NY016021"
NY016021$patient_id <- "NY016021"

NY016022T1$ICB <- "aPD1"
NY016022T1$Sex <- "F"
NY016022T1$Age <- "66"
NY016022T1$ICB_response <- "MPR"
NY016022T1$Stage <- "II"
NY016022T1$sample_id <- "NY016022T1"
NY016022T1$patient_id <- "NY016022"

NY016022T2$ICB <- "aPD1"
NY016022T2$Sex <- "F"
NY016022T2$Age <- "66"
NY016022T2$ICB_response <- "MPR"
NY016022T2$Stage <- "II"
NY016022T2$sample_id <- "NY016022T2"
NY016022T2$patient_id <- "NY016022"

NY016022T3$ICB <- "aPD1"
NY016022T3$Sex <- "F"
NY016022T3$Age <- "66"
NY016022T3$ICB_response <- "MPR"
NY016022T3$Stage <- "II"
NY016022T3$sample_id <- "NY016022T3"
NY016022T3$patient_id <- "NY016022"

NY016022T4$ICB <- "aPD1"
NY016022T4$Sex <- "F"
NY016022T4$Age <- "66"
NY016022T4$ICB_response <- "MPR"
NY016022T4$Stage <- "II"
NY016022T4$sample_id <- "NY016022T4"
NY016022T4$patient_id <- "NY016022"

NY016025T1$ICB <- "aPD1"
NY016025T1$Sex <- "F"
NY016025T1$Age <- "74"
NY016025T1$ICB_response <- "MPR"
NY016025T1$Stage <- "III"
NY016025T1$sample_id <- "NY016025T1"
NY016025T1$patient_id <- "NY016025"

NY016025T2$ICB <- "aPD1"
NY016025T2$Sex <- "F"
NY016025T2$Age <- "74"
NY016025T2$ICB_response <- "MPR"
NY016025T2$Stage <- "III"
NY016025T2$sample_id <- "NY016025T2"
NY016025T2$patient_id <- "NY016025"

NY016025T3$ICB <- "aPD1"
NY016025T3$Sex <- "F"
NY016025T3$Age <- "74"
NY016025T3$ICB_response <- "MPR"
NY016025T3$Stage <- "III"
NY016025T3$sample_id <- "NY016025T3"
NY016025T3$patient_id <- "NY016025"

NY016025T4$ICB <- "aPD1"
NY016025T4$Sex <- "F"
NY016025T4$Age <- "74"
NY016025T4$ICB_response <- "MPR"
NY016025T4$Stage <- "III"
NY016025T4$sample_id <- "NY016025T4"
NY016025T4$patient_id <- "NY016025"

#Merge objects and create combined seurat obj

NY <- merge(NY016016, y = c(NY016021, NY016007T1, NY016007T2, NY016007T3, NY016014T1, NY016014T2, NY016014T3,
                           NY016015T1, NY016015T2, NY016015T3, NY016015T4, NY016022T1, NY016022T2, NY016022T3,
                           NY016022T4, NY016025T1, NY016025T2, NY016025T3, NY016025T4), project = "NY",
            add.cell.ids = c("NY016016", "NY016021", "NY016007T1", "NY016007T2", "NY016007T3", "NY016014T1",
                             "NY016014T2", "NY016014T3", "NY016015T1", "NY016015T2", "NY016015T3", "NY016015T4",
                             "NY016022T1", "NY016022T2", "NY016022T3", "NY016022T4", "NY016025T1", "NY016025T2",
                             "NY016025T3", "NY016025T4"))
                             
#QC estimation

NY[["percent.mt"]] <- PercentageFeatureSet(NY, pattern = "^MT-")
NY[["percent.hb"]] <- PercentageFeatureSet(NY, pattern = "^HBA|^HBB")
NY[["percent.rp"]] <- PercentageFeatureSet(NY, pattern = "^RPS|^RPL")


nFeature_lower <- 500
nFeature_upper <- 10000
nCount_lower <- 1000
nCount_upper <- 100000
percent.mt_lower <- 0
percent.mt_upper <- 30
percent.hb_lower <- 0
percent.hb_upper <- 5


#Filter low quality cells
NY <- subset(NY, subset = nFeature_RNA > nFeature_lower & nFeature_RNA < nFeature_upper & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper & percent.mt < percent.mt_upper & percent.hb < percent.hb_upper)

#Prepare metadata and counts matrix for Scanpy
NY$barcode <- colnames(NY)
write.csv(NY@meta.data, file='D:Scanpy/new_york/metadata.csv', quote=F, row.names=F)

counts_matrix <- GetAssayData(NY, assay='RNA', slot='counts')
writeMM(counts_matrix, file='D:Scanpy/new_york/counts.mtx')

write.table(
  data.frame('gene'=rownames(counts_matrix)),file='D:Scanpy/new_york/genes.csv',
  quote=F,row.names=F,col.names=F
)
