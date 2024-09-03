# markers by set
Set_I <- c("mic_SMMS2", "mic_SMMS13", "mic_SMMS16","Shae10", "D105", "C114",
           "Sha104176")
Set_II<- c("mic_SMMS3", "mic_SMMS17", "mic_SMMS18", "mic_SMMS21","C102",
           "Shae12", "C112")
Set_III<- c("C146", "Shae14", "Shae05")
Set_IV<- c("mic_13TAGA", "mic_SM13_410", "mic_1F8A", "mic_SMDA23")
Set_V<- c("mic_29E6A", "mic_SM13_478", "mic_SMU31768", "mic_15J15A")
Set_VI<- c("mic_LG5_sc36b ", "mic_sc23b", "mic_SMD28")
Set_VII<- c("mic_L46951", "mic_R95529", "mic_LG1_sc276", "mic_LG5_sc475")

mic_set_list <- list(Set_I, Set_II, Set_III, Set_IV, Set_V, Set_VI, Set_VII)
names(mic_set_list) <- c("Set_I", "Set_II", "Set_III", "Set_IV", "Set_V",
                         "Set_VI", "Set_VII")

# channel 1; 6-FAM (blue)
mic_SMMS2 <- c(211, 215, 219, 223, 227, 231, 235, 239)
mic_SMMS3 <- c(177, 180, 183, 186, 189, 192, 195, 198, 201, 204, 207)
mic_13TAGA <- c(102, 106, 110, 114, 118, 122, 126, 130, 134, 138)
mic_29E6A <- c(153, 156, 159, 162, 165, 168, 171, 174,177,180,183, 186)
mic_SM13_478 <- c(224, 227, 230, 233, 236, 239, 242, 245, 248, 251, 254, 257)
mic_SMD28 <- c(225, 228, 231, 234, 237, 240)
mic_L46951 <- c(156, 159, 162, 165, 168, 171, 174, 177, 180, 183, 186, 189)
mic_LG1_sc276 <- c(98, 101, 104, 107, 110)
Shae10 <- c(161,164,167,170,173,176,179,182,185,188,191,194,
            197,200,203,206,209,212,215,218)
C112 <- c(282,285,288,291,294,297,300,303,306,309,312,315,
          318,321,324,327,330,333,336,339)
Shae05 <- c(247,250,253,256,259,262,265,268,271,274,277,280,
            283,286,289,292,295,298,301,304)
ch1_lists <- c("mic_SMMS2", "mic_SMMS3", "mic_13TAGA", "mic_29E6A",
               "mic_SM13_478", "mic_SMD28", "mic_L46951", "mic_LG1_sc276",
               "Shae10", "C112", "Shae05")

# channel 2; VIC (green)
mic_SMMS16 <- c(210, 213, 216, 219, 222, 225, 228, 231, 234)
mic_SMMS17 <- c(286, 289, 292, 295, 298, 301, 304, 307)
mic_SMMS18 <- c(195, 198, 201, 204, 207, 210, 213, 216, 219, 222, 225, 228,231)
mic_1F8A <- c(149, 152, 155, 158, 161, 164, 167, 170)
mic_15J15A <- c(208, 211, 214, 217, 220, 223, 226, 229, 232, 235, 238, 241)
D105 <- c(161,165,169,173,177,181,185,189,193,197,201,205,
          209,213,217,221,225,229,233,237)
Shae12 <- c(226,229,232,235,238,241,244,247,250,253,256,259,
            262,265,268,271,274,277,280,283)
C146 <- c(135,138,141,144,147,150,153,156,159,162,165,168,171,
          174,177,180,183,186,189,192)
ch2_lists <- c("mic_SMMS16", "mic_SMMS17", "mic_SMMS18", "mic_1F8A",
               "mic_15J15A","D105", "Shae12", "C146")

# channel 3; NED (yellow)
mic_SMDA23 <- c(187,191, 195, 199, 203, 207, 211, 215, 219, 223, 227, 231, 235)
mic_SMU31768 <- c(185, 188, 191, 194, 197, 200, 203, 206, 209, 212, 215, 218,
                  221, 224, 227)
mic_LG5_sc36b <- c(232, 235, 238, 241, 244, 247, 250, 253, 256, 259, 262, 265,
                   268)
mic_sc23b <- c(191, 194, 197, 200, 203, 206, 209, 212, 215, 218)
mic_R95529 <- c(219, 222, 225, 228, 231, 234, 237, 240, 243, 246, 249, 252)
C114 <- c(208,211,214,217,220,223,226,229,232,235,238,241,244,
          247,250,253,256,259,262,265, 268,271,274,277,280)
Shae14 <- c(170,174,178,182,186,190,194,198,202,206,210,214,
            218,222,226,230,234,238,242,246)
ch3_lists <-c("mic_SMDA23", "mic_SMU31768", "mic_LG5_sc36b", "mic_sc23b",
              "mic_R95529","C114", "Shae14")

# channel 4; PET (red)
mic_SMMS13 <- c(183, 186, 189, 192, 195, 198, 201, 204)
mic_SMMS21 <- c(172, 175, 178, 181, 184, 187)
mic_SM13_410 <- c(191, 194, 197, 200, 203)
mic_LG5_sc475 <- c(281, 284, 287, 290, 293, 296, 299, 302, 305)
Sha104176 <- c(283,286,289,292,295,298,301,304,307,310,313,316,
               319,322,325,328,331,334,337,340)
C102 <- c(164,167,170,173,176,179,182,185,188,191,194,197,200,
          203,206,209,212,215,218,221)
ch4_lists <- c("mic_SMMS13", "mic_SMMS21", "mic_SM13_410", "mic_LG5_sc475",
               "Sha104176", "C102")

dye_mic_list <- list(ch1_lists, ch2_lists, ch3_lists, ch4_lists)
names(dye_mic_list) <- c("ch1_lists", "ch2_lists", "ch3_lists", "ch4_lists")


test_that("score_markers_rev3 returns expected output", {

  file_path <- system.file("extdata", package = "pooledpeaks")

  mock_fsa_batch_imp_output<- fsa_batch_imp(file_path, channels = 5,
                                            fourier = TRUE, saturated = TRUE,
                                            lets.pullup = FALSE,
                                            plotting = FALSE, rawPlot = FALSE,
                                            llength = 3000, ulength = 80000)
  names(mock_fsa_batch_imp_output)<-c("23.2a_I_A01_2012-07-18.fsa",
                                      "23.2b_I_A07_2012-07-18.fsa",
                                      "30.3a_I_B01_2012-07-18.fsa",
                                      "30.3b_I_B07_2012-07-18.fsa",
                                      "33.1a_I_C01_2012-07-18.fsa",
                                      "33.1b_I_C07_2012-07-18.fsa",
                       "Multiplex_set_I_Shaem.1a_1_Sample_20221028_215632.fsa",
                       "Multiplex_set_I_Shaem.1b_1_Sample_20221028_232301.fsa",
                       "Multiplex_set_I_Shaem.3a_2_Sample_20221028_215633.fsa",
                       "Multiplex_set_I_Shaem.3b_2_Sample_20221028_232302.fsa",
                       "Multiplex_set_I_Shaem.4a_3_Sample_20221028_215634.fsa",
                       "Multiplex_set_I_Shaem.4b_3_Sample_20221028_232303.fsa")

  # Dummy data and parameters
  panel <- Shae10
  ladder <- c(20, 40, 60, 80, 100, 114, 120, 140, 160, 180, 200, 214, 220,
              240, 250, 260, 280, 300, 314, 320, 340, 360, 380, 400, 414,
              420, 440, 460, 480, 500, 514, 520, 540, 560, 580, 600)

  mock_fsa_batch_imp_output <- associate_dyes(mock_fsa_batch_imp_output,
                                              file_path)

  Fragman::ladder.info.attach(stored = mock_fsa_batch_imp_output,
                              ladder = ladder,
                              ladd.init.thresh = 200, prog = FALSE,
                              draw = FALSE)

  # Ensure panel is numeric
  panel <- as.numeric(panel)

  # Running the score_markers_rev3 function with a smaller subset
  small_subset <- mock_fsa_batch_imp_output[1:6]
  # Using a smaller subset for testing

  result <- score_markers_rev3(my.inds = small_subset,
                               channel = 1,
                               channel.ladder = 5,
                               panel = "panel",
                               ladder = ladder,
                               init.thresh = 100,
                               ploidy = length(panel),
                               shift = 1,
                               windowL = 1,
                               windowR = 0.5,
                               left.cond = c(0, 2.5),
                               right.cond = 0,
                               pref = 1,
                               plotting = FALSE)

  expect_type(result, "list")
  expect_equal(length(result), 6)
  expect_true(all(sapply(result, function(x) all(c("pos", "hei", "wei") %in%
                                                   names(x)))))
})

test_that("score_markers_rev3 handles plotting option", {

  file_path <- system.file("extdata", package = "pooledpeaks")

  mock_fsa_batch_imp_output<- fsa_batch_imp(file_path, channels = 5,
                                            fourier = TRUE, saturated = TRUE,
                                            lets.pullup = FALSE,
                                            plotting = FALSE, rawPlot = FALSE,
                                            llength = 3000, ulength = 80000)
  names(mock_fsa_batch_imp_output)<-c("23.2a_I_A01_2012-07-18.fsa",
                                      "23.2b_I_A07_2012-07-18.fsa",
                                      "30.3a_I_B01_2012-07-18.fsa",
                                      "30.3b_I_B07_2012-07-18.fsa",
                                      "33.1a_I_C01_2012-07-18.fsa",
                                      "33.1b_I_C07_2012-07-18.fsa",
                       "Multiplex_set_I_Shaem.1a_1_Sample_20221028_215632.fsa",
                       "Multiplex_set_I_Shaem.1b_1_Sample_20221028_232301.fsa",
                       "Multiplex_set_I_Shaem.3a_2_Sample_20221028_215633.fsa",
                       "Multiplex_set_I_Shaem.3b_2_Sample_20221028_232302.fsa",
                       "Multiplex_set_I_Shaem.4a_3_Sample_20221028_215634.fsa",
                       "Multiplex_set_I_Shaem.4b_3_Sample_20221028_232303.fsa")

  # Dummy data and parameters
  panel <- Shae10
  ladder <- c(20, 40, 60, 80, 100, 114, 120, 140, 160, 180, 200, 214, 220,
              240, 250, 260, 280, 300, 314, 320, 340, 360, 380, 400, 414,
              420, 440, 460, 480, 500, 514, 520, 540, 560, 580, 600)

  mock_fsa_batch_imp_output <- associate_dyes(mock_fsa_batch_imp_output,
                                              file_path)

  Fragman::ladder.info.attach(stored = mock_fsa_batch_imp_output,
                              ladder = ladder,
                              ladd.init.thresh = 200, prog = FALSE,
                              draw = FALSE)

  # Ensure panel is numeric
  panel <- as.numeric(panel)

  # Test without plotting
  result <- score_markers_rev3(my.inds = mock_fsa_batch_imp_output,
                               channel = 1,
                               channel.ladder = 5,
                               panel = "panel",
                               ladder = ladder,
                               init.thresh = 100,
                               ploidy = length(panel),
                               shift = 1,
                               windowL = 1,
                               windowR = 0.5,
                               left.cond = c(0, 2.5),
                               right.cond = 0,
                               pref = 1,
                               plotting = FALSE)
  expect_type(result, "list")

  # Test with plotting
  plot_dir <- file.path(tempdir(), "plots_scoring")
  if (Sys.info()["sysname"] == "Windows") {
    # On Windows, use backslash as the separator
    plot_dir <- paste0(plot_dir, "\\")
  } else {
    # On Unix-like systems, use forward slash
    plot_dir <- paste0(plot_dir, "/")
  }

  result <- score_markers_rev3(my.inds = mock_fsa_batch_imp_output,
                               channel = 1,
                               channel.ladder = 5,
                               panel = "panel",
                               ladder = ladder,
                               init.thresh = 100,
                               ploidy = length(panel),
                               shift = 1,
                               windowL = 1,
                               windowR = 0.5,
                               left.cond = c(0, 2.5),
                               right.cond = 0,
                               pref = 1,
                               plotting = TRUE,
                               plotdir = plot_dir)
  expect_type(result, "list")
  expect_true(file.exists(file.path(plot_dir,"all_panel_scores.pdf")))
  unlink(plot_dir, recursive = TRUE)
})
