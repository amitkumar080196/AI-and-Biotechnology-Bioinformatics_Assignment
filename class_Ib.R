dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")
data <- read.csv(file.choose())
View(data)
str(data)

disease_status <- c("cancer","normal", "cancer", "normal")
class(disease_status)
disease_status <- as.factor(disease_status)
class(disease_status)
disease_status


disease_status <- factor(disease_status,
                         levels = c("normal", "cancer"))
disease_status

data$gender_fac <- as.factor(data$gender)
str(data)

data$gender_num <- ifelse(data$gender_fac == "female", 1, 0)
data$gender_num <- as.factor(data$gender_num)

class(data$gender_num)

data$smoker_fac <- as.factor(data$smoker)
str(data)

data$smoker_num <- ifelse(data$smoker_fac == "yes", 2, 0)
data$smoker_num <- as.factor(data$smoker_num)


mean(bmi)
mean_results <- mean(bmi)
plot(bmi)
hist(bmi)
barplot(bmi)
summary(bmi)






write.csv(disease_status, file = "results_1/patient_data.csv")

save.image (file = "full_workspace.RData")






