library(readstata13)
data_dir <- paste0("~/Dropbox/GPH_ExaminerDesign", "/Applications/STAR/Data/STARgk_Lambdas.dta")
dt <- readstata13::read.dta13(data_dir, generate.factors=TRUE)
dt$treatment <- "regular"
dt$treatment[dt$gksmall==1] <- "small"
dt$treatment[dt$gkaide==1] <- "aide"
dt$treatment <- factor(dt$treatment, levels=c("regular", "small", "aide"))
dt <- dt[, c("gkaggregate2", "treatment", "gkschid", "gktchid")]
names(dt) <- c("score", "treatment", "school", "teacher")
dt$school <- as.factor(dt$school)
dt$teacher <- as.factor(dt$teacher)

Wm <- outer(dt$school, levels(dt$school), `==`) + 0
Y  <- dt$score
X  <- dt$treatment
rg <- te_estimates(dt$score, dt$treatment, Wm)
print(knitr::kable(rg))
