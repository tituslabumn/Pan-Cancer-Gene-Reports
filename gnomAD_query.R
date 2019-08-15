# gnomAD_quer.R

cat("#####################################################################\n")
cat("############# Retrieving all gnomAD variants for GOI ################\n")
cat("#####################################################################\n\n\n")

cat("Cheching for linked pcgr_selenium_chrome docker container...\n\n")

if("selenium-file" %in% list.files(path = "/")){
  cat("\tlink confirmed :]\n\n") #assumes selenium container was not stopped after pcgr container was brought up
}else{
  stop("\tCan't fined shared volume from selenium container")
}

cat("Retriving gnomAD non-cancer variant file via RSelenium web scraper\n\n")

  # docker run --name pcgr_selenium_chrome --rm -d -p 4444:4444 -v /home/seluser/Downloads selenium/standalone-chrome
# docker run -ti --rm --volumes-from pcgr_selenium_chrome -v ~/PCGR_OUTPUT:/OUTPUT --link=pcgr_selenium_chrome tsharding/pcgr_v1.0
# docker stop pcgr_selenium_chrome

library("RSelenium")

URL_gnomAD <- paste0("https://gnomad.broadinstitute.org/gene/",GOI_ENSG,"?dataset=gnomad_r2_1_non_cancer")
# https://gnomad.broadinstitute.org/gene/ENSG00000145555?dataset=gnomad_r2_1_non_cancer

# eCaps <- list(
#   chromeOptions = 
#     list(prefs = list(
#       "profile.default_content_settings.popups" = 0L,
#       "download.prompt_for_download" = FALSE,
#       "download.default_directory" = "/OUTPUT"
#     )
#     )
# )

# selenium is the service name of other container on user defined network per docker-compose file
remDr <- remoteDriver(remoteServerAddr = "selenium" ,port = 4444L, browser = "chrome", )#extraCapabilities = eCaps) 
cat("\tOpening connection via remote driver\n")
remDr$open(silent = TRUE)
# navigate to GOI gnomAD page
cat("\tNavigating to gnomAD browser GOI page(non-cancer):\n")
cat("\t\t",URL_gnomAD,"\n")
remDr$navigate(URL_gnomAD)
Sys.sleep(3)
cat("\t\tPage title: ",remDr$getTitle(),"\n")
      # remDr$screenshot(display = TRUE)
#   find download button (xpath found by exploring element in chromium)
cat("\tFinding download button\n")
webElem_dl_button <- remDr$findElement(using = "xpath","//*[@id='root']/div/div/div[2]/div/div[5]/section/div[2]/button")
Sys.sleep(3)
# click download button
cat("\tMouse over download button\n")
webElem_dl_button$sendKeysToElement(list(key = "down_arrow"))
# remDr$screenshot(display = TRUE)
Sys.sleep(3)
cat("\tClick download button\n")
webElem_dl_button$clickElement()
cat("\tClose driver connection\n")
remDr$close()

print(list.files(path = "/selenium-file/"))
