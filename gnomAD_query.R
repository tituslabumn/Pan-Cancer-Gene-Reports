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

remDr <- remoteDriver(remoteServerAddr = "selenium_chrome" ,port = 4444L, browser = "chrome", )#extraCapabilities = eCaps)
remDr$open(silent = TRUE)
# navigate to GOI gnomAD page
remDr$navigate(URL_gnomAD)
Sys.sleep(3)
      # remDr$getTitle()
      # remDr$screenshot(display = TRUE)
# find download button (xpath found by exploring element in chromium)
webElem_dl_button <- remDr$findElement(using = "xpath","//*[@id='root']/div/div/div[2]/div/div[5]/section/div[2]/button")
Sys.sleep(3)
# click download button
webElem_dl_button$sendKeysToElement(list(key = "down_arrow"))
# remDr$screenshot(display = TRUE)
Sys.sleep(3)
webElem_dl_button$clickElement()
remDr$close()

print(list.files(path = "/selenium-file/"))
