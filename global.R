# Brandon Ng
# 
# Global fields for registering files into database
# 



# Global fields  ----
proteins <- c("ERK", "PAD4", "IRAK3", "IRAK4") %>% sort()
cells <- c("THP-1", "EL4", "RA1", "CT26", "A375") %>% sort()
scientists <- c("Brandon", "Dan", "Leslie", "Gramoz", "Tom", "Vanessa") %>% sort()
drugs <- c("123","456","12424") %>% sort()
fields <- c("form_data","date", "scientist", "protein", "cell", "description","drug")
empty_fields <- c("form_data","scientist", "protein", "cell", "description","drug")