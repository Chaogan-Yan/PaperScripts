# 统计森林图
library("metansue")
library("R.matlab")
library("stringr")
library("ggplot2")
library(R.matlab)

subs <- c("AdultSubjects")
for (sub in subs) {
    LR_all <- c("LH","RH")
    model_all <- c("")
    path_data_all <- paste0('',sub,'')

  
    for (LR in LR_all){
        n_num <- 1
        for (model in model_all){
            N <- readMat(paste0(path_data_all,model,'_ForR.mat')) #用于后面的meta分析的N值
            path_all <- '' #提取出的每个站点的值的路径
            value_name <- paste0(path_all, 'peak_value/',sub,'/',model,'_',LR,'extract.csv')
            if(!file.exists(value_name)) {
                cat(value_name, "does not exist. Skipping...\n")
                next
            }
            value <- read.csv(value_name,header=F) #提取的值
            save_path <- paste0(path_all,'forest_plot/',sub,'/',model,'/') #森林图存放
            save_path_new <- paste0(path_all,'forest_plot_new/',sub,'/',model,'/') # 加了标题的森林图
            result_save_path <- paste0(path_all, 'result/',sub,'/') #存放统计量结果

            if (!dir.exists(save_path)) {
                # 文件夹不存在，创建它
                dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
                cat("文件夹已创建成功！\n")
                } else {
                cat("文件夹已存在！\n")
                }
            if (!dir.exists(save_path_new)) {
                # 文件夹不存在，创建它
                dir.create(save_path_new, recursive = TRUE, showWarnings = FALSE)
                cat("文件夹已创建成功！\n")
                } else {
                cat("文件夹已存在！\n")
                }
            if (!dir.exists(result_save_path)) {
                # 文件夹不存在，创建它
                dir.create(result_save_path, recursive = TRUE, showWarnings = FALSE)
                cat("文件夹已创建成功！\n")
                } else {
                cat("文件夹已存在！\n")
                }
            roi_num <- ncol(value)-1
            if (roi_num == 0) {
                cat("roi_num is 0. Skipping...\n")
                next
            }
            site_num <- nrow(value)-2
            ROI_list <- value[2,1:roi_num+1]
            site_list <- value[1:site_num+2,1]
            


            # 初始化指标
            n = 1
            Z = matrix(0,roi_num, 1)
            Cohend = matrix(0,roi_num, 1)
            se = matrix(0,roi_num, 1)
            ci = matrix("",roi_num, 1)
            p = matrix(0,roi_num, 1)
            i2 = matrix(0,roi_num, 1)
            ROI = matrix("",roi_num, 1)

            #对每个结果图
            for (i in 1:roi_num){
                N1 <- N$N1[,1]
                if(model == "M9_HDRS" || model == "M2_DxByAge"){
                    N2 <- 0
                }else{
                    N2 <- N$N2[,1]
                }
                # 读取结果图
                if(length(N2) == 1){
                    x <- smc_from_t(c(as.numeric(value[1:site_num+2, i+1])),c(N1),labels = site_list)
                } else {
                    x <- smd_from_t(c(as.numeric(value[1:site_num+2, i+1])),c(N1),c(N2),labels = site_list)
                }
                m <- meta(x)
            
                fig <- forest(m)

                output_string <- capture.output(m)
                words <- strsplit(output_string[20], " ") # 捕获m的输出，并用于后面的CI获取

                Cohend[n] = m[['hypothesis']][['coef']]
                Z[n] = m[['hypothesis']][['z']]
                se[n] = m[['hypothesis']][['se']]
                if (words[[1]][14] != ""){
                    ci[n] = paste0("[",words[[1]][12],", ",words[[1]][14],"]")
                } else if(words[[1]][15] != ""){
                    ci[n] = paste0("[",words[[1]][12],", ",words[[1]][15],"]")
                } else if(words[[1]][13] != ""){
                    ci[n] = paste0("[",words[[1]][12],", ",words[[1]][13],"]")
                }
                # ci[n] = paste0("[",words[[1]][12],", ",words[[1]][14],"]")
                p[n] = m[['hypothesis']][['p.value']]
                i2[n] =  m[['heterogeneity']][['i2']]
                if (length(ROI_list) == 1){
                    ROI[n] = paste0(ROI_list)
                } else{
                    ROI[n] = paste0(ROI_list[1, i])
                }
                fig_name <- paste0(model, "_", LR, "_", i)
                # 创建一个PNG文件
                png(filename = paste0(save_path,fig_name,".png"), width = 1200, height = 1800, res = 130)
                par(mar = c(1, 1, 1, 1))
                # 画图
                print(forest(m))
                # 结束并保存PNG文件
                dev.off()
                n = n + 1
            }
            df <- data.frame(
                "ROI" = ROI,
                "Side" = LR,
                "Z" = Z,
                "Cohen's d " = Cohend,
                "Std.Err" = se,
                "CI" = ci,
                #   "% Difference" = 
                "p-value" = p,
                "I2" =   i2
                )
                # tab_name <- paste(split_str[1:(length(split_str)-3)], collapse = "-")
                tab_name <- paste0(model, "_", LR, "_stat")

                # 保存数据框为 CSV 文件
                write.table(df, file = paste0(result_save_path,tab_name,".csv"), sep = ",", row.names = FALSE)
                # 保存数据框为 RData 文件

        }
    }
}
