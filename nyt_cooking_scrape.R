


library(tidyverse)
library(httr)
library(rvest)

root_url <- 'https://cooking.nytimes.com/'

nyt_cooking <- read_html(root_url)

root_links <- nyt_cooking %>% html_nodes('.recipe-card') %>% 
  html_nodes("div") %>%
  html_nodes("a") %>%
  html_attr(., 'href') %>% unique

root_links <- root_links[str_detect(root_links, "recipes/")]

other_page <- read_html(paste0(root_url, root_links[2]))

page_tags <- other_page %>% html_nodes(".recipe-metadata") %>% html_nodes(".tags-nutrition-container") %>%
  html_nodes(".tag") %>% html_text()

other_page %>% html_nodes(".recipe-metadata") %>% html_nodes(".recipe-actions") %>%
  html_nodes("#recipe-ratings") 



new_links <- other_page %>% html_nodes('.recipe-card') %>% 
  html_nodes("div") %>%
  html_nodes("a") %>%
  html_attr(., 'href') %>% unique







