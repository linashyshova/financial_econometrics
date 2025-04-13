library(ggplot2)
library(scales)

# Read data
df <- read.table("assignment_1/data/HSI.txt")

# Add column name
names(df) <- "price"    

# Part A.1

# Add another column to indicate the date. Let's assume we started collecting data from 01.01.1997
start_date <- as.Date("1997-01-01")
df$date <- seq(from = start_date, by = "week", length.out = nrow(df))

# Plot weekly price
ggplot(df, aes(x = date, y = price)) +
  geom_line(color = "blue") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = "Year", y = "Price USD", title = "HSI Index Price - Weekly") +
  theme(axis.text.x = element_text(angle = 90))

df$log_return <- c(NA, 100 * diff(log(df$price)))

# Plot weekly log returns
ggplot(df, aes(x = date, y = log_return)) +
  geom_line(na.rm = TRUE, color = "blue") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = "Year", y = "Log Return (%)", title = "Weekly Log Returns of Stock Price") +
  theme(axis.text.x = element_text(angle = 90))

# ACF for weekly log returns
acf(df$log_return, na.action = na.pass, main = "")

# ACF for weekly squared log returns
acf(df$log_return^2, na.action = na.pass, main = "")
