if (!require(httr)) install.packages("httr")
if (!require(desc)) install.packages("desc")
library(httr)
library(desc)
url <- "https://raw.githubusercontent.com/spang-lab/FastRet/main/DESCRIPTION"
local_file <- "./DESCRIPTION"
response <- GET(url)
if (http_status(response)$category != "Success") stop("Failed to get DESCRIPTION file from URL")
desc_from_url <- desc::desc(text = content(response, "text"))
desc_local <- desc::desc(local_file)
version_from_url <- desc_from_url$get_version()
version_local <- desc_local$get_version()
cat("Version from URL: ", as.character(version_from_url), "\n")
cat("Version from local: ", as.character(version_local), "\n")
if (version_local > version_from_url) {
    cat("Version has been increased ==> check passed\n")
    quit(status = 0, save = "no")
} else {
    stop("Version has not been increased ==> check failed\n")
    quit(status = 1, save = "no")
}
