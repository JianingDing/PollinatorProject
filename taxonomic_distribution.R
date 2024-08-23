# load required libraries
library(raster)
library(ggplot2)
library(dplyr)
library(data.table)
library(rworldmap) 
library(rworldxtra)
library(lme4)
library(cowplot)
library(patchwork)
library(viridis)
library(knitr)
library(raster)
library(ggplot2)
library(dplyr)
library(data.table)
library(rworldmap) 
library(rworldxtra)
library(lme4)
library(cowplot)
library(viridis)
library(ncdf4)
library(nlme)
library(sf)
# read in the original predicts database 
PREDICTS <- readRDS("C:\\Users\\dingj\\Desktop\\project\\predict data_use\\PREDICTS_pollinators_8_exp.rds") %>%
  mutate(COL_ID = as.character(COL_ID))

# source in additional functions
source("00_functions.R")
source("CorrectSamplingEffort.R")
source("SiteMetrics.R")
source("MergeSites.R")
source("ReadPREDICTS.R")

# calc baseline (mean and sd)
calc_baseline <- function(data_file, func, pred_points, pred_points_sp){
  
  # calcualte either the mean or standard error for baseline, then extract points for predicts sites
  data_fin <- calc(data_file, func) %>%
    extract(pred_points_sp)
  
  # bind the extracted values back onto the predicts coordinates
  data_fin <- data.frame(pred_points[,1:3 ], data_fin)
  
  return(data_fin)
  
}

# Load in the mean temperature data from CRU
tmp <- raster::stack("C:\\Users\\dingj\\Desktop\\project\\climate data\\cru_ts4.03.1901.2018.tmp.dat.nc", varname="tmp")

# Read in the predicts pollinators
PREDICTS_pollinators_orig <- readRDS("C:\\Users\\dingj\\Desktop\\project\\predict data_use\\PREDICTS_pollinators_8_exp.rds") %>%
  dplyr::select(-clade_rank, -confidence)  %>%
  filter(Class == "Insecta") %>%
  dplyr::filter(Predominant_land_use %in% c("Cropland", "Primary vegetation")) %>%
  droplevels()

# Create table of number of unique species
PREDICTS_pollinators_orig %>%
  select(Class, Order, Best_guess_binomial) %>%
  filter(Best_guess_binomial != "") %>%
  unique() %>%
  group_by(Class, Order) %>%
  tally()

# Number of unique sites per land use
PREDICTS_pollinators_orig %>%
  select(Predominant_land_use, SSBS) %>%
  unique() %>%
  group_by(Predominant_land_use) %>%
  tally()

# Calculate variation in sampling effect
PREDICTS_pollinators_orig %>%
  group_by(SS) %>%
  summarise(sampling_variation = sd(Sampling_effort, na.rm = TRUE)) %>%
  mutate(varies = ifelse(sampling_variation > 0, 1, 0)) %>%
  mutate(rows = 1) %>%
  summarise(total_vary = sum(varies, na.rm = TRUE), total = sum(rows)) %>%
  mutate(percentage = (total_vary / total) * 100)

# Set up vector for filtering for vertebrates and invertebrates
predict_climate_list <- list()
PREDICTS_pollinators_orig <- list(PREDICTS_pollinators_orig)

# PREDICTS data compilation
# filter for main pollinating taxa
PREDICTS_pollinators <- PREDICTS_pollinators_orig[[1]]

# correct for sampling effort 标准化采样努力，以便不同站点或记录之间可以公平比较
PREDICTS_pollinators <- predictsFunctions::CorrectSamplingEffort(PREDICTS_pollinators)

# calculate site metrics including all species (confirmed and not confirmed pollinator)
order.sites.div <-predictsFunctions::SiteMetrics(diversity = PREDICTS_pollinators,
                                                 extra.cols = c("SSB", "SSBS", "Predominant_land_use", "UN_region", "Pollinating"))


# set id column for merging back into correct place设置合并用的ID列
order.sites.div$id_col <- 1:nrow(order.sites.div)  

# PREDICTS sites with the month of the recording提取和转换日期
PRED_sites <- order.sites.div %>% select(id_col, Latitude, Longitude, Sample_end_latest) %>%
  mutate(Sample_end_latest = paste("X", substr(Sample_end_latest, start = 1, stop = 7), sep = "")) %>%
  mutate(Sample_end_latest = gsub("-", ".", Sample_end_latest)) %>%
  filter(!is.na(Latitude))  #去除那些纬度数据为缺失的记录

# calculate the means and standard deviation for the beginning of the series
# take names of values for 1901 to 1930
tmp1901_1930 <- tmp[[names(tmp)[1:360]]]  

# extract the points for each the predicts coordinates提取地理坐标并创建空间点
PRED_sites_sp <- PRED_sites %>%
  select(Longitude, Latitude) %>%
  filter(!is.na(Latitude)) %>%
  SpatialPoints() #用于根据经纬度数据创建空间点对象

# calculate the mean baseline, and convert to character for merging 计算平均基线值并转换为字符类型
climate_start_mean <- calc_baseline(tmp1901_1930, 
                                    func = mean, 
                                    pred_points = PRED_sites, #站点数据
                                    pred_points_sp = PRED_sites_sp) %>% #空间点数据
  mutate(Latitude = as.character(Latitude)) %>%
  mutate(Longitude = as.character(Longitude))

# calculate the sd baseline, and convert to character for merging 计算标准偏差基线值并转换为字符类型
climate_start_sd <- calc_baseline(tmp1901_1930, 
                                  func = stats::sd, 
                                  pred_points = PRED_sites, 
                                  pred_points_sp = PRED_sites_sp) %>%
  mutate(Latitude = as.character(Latitude)) %>%
  mutate(Longitude = as.character(Longitude))

#以上代码：处理和分析基于地理位置的气候数据，计算出特定时间段内的平均气候基线和气候变异性

# calculate the mean temperatures for each predicts site, 11 months previously
# set up empty list for each dataframe
raster_means <- list()  

# time the length of the loop over each unique date in predicts sites
system.time(
  # for each unique site month, select the pixels for that date, and 11 months previously, and then convert each raster set to a dataframe
  # 计算每个预测站点在过去11个月内的平均温度
  
  for(i in 1:length(unique(PRED_sites$Sample_end_latest))){
    
    # create unique list of the end dates for the predicts sites
    pred_dates <- unique(PRED_sites$Sample_end_latest)
    
    # select the raster pixels for the end date of the predicts site, and then the index of that name
    date_raster <- grepl(pred_dates[i], names(tmp))#检查 tmp 数据集（可能是一个包含多个栅格层的列表或数组）的名称中是否包含 pred_dates[i]（当前循环的日期）。grepl 返回一个逻辑向量，表明哪些元素匹配。
    site_index <- which(date_raster == TRUE)#找出上述逻辑向量中为 TRUE 的索引，这些索引指向与当前日期匹配的栅格层。
    
    # select the previous 11 indices - i.e. 11 months previous worth of pixels
    ind_raster <- tmp[[names(tmp)[(site_index - 11): site_index]]] #使用这些索引来从 tmp 中选择对应的栅格数据，这包括了从当前月份往前数的11个月的数据。names(tmp)[(site_index - 11): site_index] 计算得到需要的栅格层的名称，然后从 tmp 中提取这些层。
    
    # filter the raster for only coordinates we have PREDICTS sites for that date
    # 筛选特定日期的预测站点并创建空间点
    PRED_sites_filt <- PRED_sites %>%
      filter(Sample_end_latest == pred_dates[i]) %>%
      select(Longitude, Latitude) %>%
      SpatialPoints() #将筛选后的经度和纬度转换为空间点对象，这对于后续的空间数据分析非常重要。
    
    # identify site ids for merging 识别站点ID以便合并
    site_ids <- PRED_sites %>%
      filter(Sample_end_latest == pred_dates[i]) %>%
      select(id_col, Longitude, Latitude)
    
    # filter the raster for that date for the locations we have predicts sites根据预测站点位置筛选栅格数据
    PRED_coords <- cbind(site_ids, extract(ind_raster, PRED_sites_filt, na.rm = FALSE))
    
    # convert that set of dates to a dataframe
    ind_raster_frame <- as.data.frame(PRED_coords)
    
    # remove the extra coordinate columns for calculating the row means 移除额外的坐标列
    ind_raster_values <- ind_raster_frame %>% select(-id_col, -Longitude, -Latitude)
    
    # calculate the mean values for each coordinate and bind back onto the coordinates 计算每个坐标的平均值并重新绑定
    raster_means[[i]] <- cbind(names(tmp)[site_index], (ind_raster_frame %>% select(id_col, Longitude, Latitude)), rowMeans(ind_raster_values))
    colnames(raster_means[[i]]) <- c("end_date", "id_col", "x", "y", "mean_value")
    
    # print the iteration number
    print(i)
  }
  #目的是为每个独特的样本结束日期计算相关预测站点的平均温度，并将这些信息与站点的地理坐标一起存储
)  

## adjust the mean value for each site for the baseline at that site
#如何调整每个站点的平均温度值，使其与该站点的基线（即早期记录的平均值和标准偏差）相对应。
#这是一个重要的分析步骤，用于计算温度异常和标准化温度异常
#可以帮助了解气候变化的具体影响。
# first, merge the baseline sd and mean by coordinate for each site
adjusted_climate <- rbindlist(raster_means) %>%
  select(-end_date) %>%
  unique() %>%
  inner_join(climate_start_mean, by = "id_col") %>%#为每个站点添加了基线平均温度 
  rename("mean_base" = "data_fin") %>%
  inner_join(climate_start_sd, by = "id_col") %>% #再次进行内连接，这次是添加基线标准偏差 
  rename("sd_base" = "data_fin") %>%
  mutate(anomaly = mean_value - mean_base) %>% #表示温度相对于历史平均的偏差
  mutate(standard_anom = anomaly / sd_base) #标准化的异常值表示每个站点的温度偏差与历史变异性的比例，有助于在不同站点间进行比较

# bind the adjusted climate data back onto the predicts sites 
# 将气候异常值和标准化异常值与其他生态和地理数据结合，为后续的生态分析提供一个综合的数据集
predicts_climate <- inner_join(order.sites.div, adjusted_climate, by = "id_col")


# assign to list of predicts_climate and insects and vertebrates
predict_climate_list[[1]] <- predicts_climate



# build map
base_map <- get_basemap()

# fortify the main map
map_fort <- fortify(base_map)

# build map for distribution of sites
site_distribution_pollinator <- ggplot() + 
  geom_polygon(aes(x = long, y = lat, group = group), data = map_fort, fill = "lightgrey") +
  geom_point(aes(x = Longitude, y = Latitude, fill = Predominant_land_use), data = predict_climate_list[[1]], shape = 21, colour = "black", alpha = 0.3) +
  scale_fill_manual("Predominant land-use", values = c("#009E73", "#E69F00")) +
  coord_map(projection = "mollweide") +
  labs(subtitle = "(A)") +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 13),
        panel.background = element_rect(fill = "white"),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        legend.key=element_blank(), legend.position = "bottom", 
        legend.justification = "left")

print(site_distribution_pollinator)
ggsave("site_distribution_landuse.png", plot = site_distribution_pollinator, width = 10, height = 7, dpi = 300)

# build map for distribution of sites
anomaly_pollinator <- ggplot() + 
  geom_polygon(aes(x = long, y = lat, group = group), data = map_fort, fill = "lightgrey") +
  geom_point(aes(x = Longitude, y = Latitude, colour = standard_anom), data = predict_climate_list[[1]], alpha = 0.3) +
  scale_colour_viridis("Std. temp.\nanomaly", option = "B") +
  coord_map(projection = "mollweide") +
  labs(subtitle = "(B)") +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 13),
        panel.background = element_rect(fill = "white"),
        legend.title = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        legend.key=element_blank(), legend.position = "bottom", 
        legend.justification = "left")

print(anomaly_pollinator)
ggsave("site_distribution_anomaly.png", plot = anomaly_pollinator, width = 10, height = 7, dpi = 300)

