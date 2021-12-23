# Importing IPUMS data
ddi <- read_ipums_ddi("~/Downloads/usa_00014.xml")
ACS_2011_2015 <- read_ipums_micro(ddi)

# Filtering out group quarters and extra people in each HH, selecting vars
ACS_2011_2015 <- ACS_2011_2015 %>% filter(!(GQ == 4 | GQ == 3)) %>% filter(PERNUM == 1) %>% select(PUMA, HHINCOME, HHWT)
attributes(ACS_2011_2015$PUMA) <- NULL
attributes(ACS_2011_2015$HHINCOME) <- NULL
# renaming vars
ACS_2011_2015 <- ACS_2011_2015 %>% rename(income = HHINCOME, wt = HHWT, PUMID = PUMA) 
# Changing class of income attribute to numeric (no labels)
attr(ACS_2011_2015$income, 'class') <- 'numeric'
# Changing negative and 0 incomes to $1
ACS_2011_2015 <- ACS_2011_2015 %>% mutate(income = ifelse(income <= 0, 1, income))

# Storing data as pums_data
pums_data <- ACS_2011_2015
rm(ACS_2011_2015, ddi)

# Writing to discc
write_csv(pums_data, "~/lorenz_revision3/code/pums_data.csv")
