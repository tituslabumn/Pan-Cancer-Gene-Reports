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

eCaps <- list(
  chromeOptions =
    list(prefs = list(
      "profile.default_content_settings.popups" = 0L,
      "download.prompt_for_download" = FALSE#,
      #"download.default_directory" = "/home/seluser/Downloads"
    )
    )
)

save.image("troubleshooting_workspace.RData") #####################

# selenium is the service name of other container on user defined network per docker-compose file
remDr <- remoteDriver(remoteServerAddr = "selenium" ,port = 4444L, browser = "chrome", extraCapabilities = eCaps) 
cat("\tOpening connection via remote driver\n")
remDr$open(silent = TRUE)
# navigate to GOI gnomAD page
cat("\tNavigating to gnomAD browser GOI page(non-cancer):\n")
cat("\t\t",URL_gnomAD,"\n")
remDr$navigate(URL_gnomAD)
Sys.sleep(3) #wait for page to load
cat("\t\tPage title:\n")
print(remDr$getTitle())
      # remDr$screenshot(display = TRUE)
#   find download button (xpath found by exploring element in chromium)
cat("\tFinding download button\n")
webElem_dl_button <- remDr$findElement(using = "xpath","//*[@id='root']/div/div/div[2]/div/div[5]/section/div[2]/button")
# click download button
Sys.sleep(3)
cat("\tMouse over download button\n")
webElem_dl_button$sendKeysToElement(list(key = "down_arrow"))
cat("\tClick download button\n")
Sys.sleep(3)
webElem_dl_button$clickElement()
cat("\tClose driver connection\n")
Sys.sleep(3)
remDr$close()

#check if downloads dir was created and contains one file only
if("Downloads" %in% list.files(path = "/selenium-file/") & length(list.files(path = "/selenium-file/Downloads/")) == 1){
  cat("File retrieved:")
  print(list.files(path = "/selenium-file/Downloads/"))
  gnomeAD_filename <- list.files(path = "/selenium-file/Downloads/")[1]
}else{
  stop("gnomAD file not retrived successfully")
}

save.image("troubleshooting_workspace.RData") #####################

cat("parsing in file\n")
GOI_gnomAD_df <- read.csv(paste0("/selenium-file/Downloads/",gnomeAD_filename), stringsAsFactors = FALSE)

save.image("troubleshooting_workspace.RData") #####################

cat("removing file dir\n")
unlink("/selenium-file/Downloads/", recursive = TRUE)
print(str(GOI_gnomAD_df))

# this may already be filtered
cat("\tFiltering out variants that lack 'PASS' annotation\n")
GOI_gnomAD_df_filtered <- GOI_gnomAD_df[GOI_gnomAD_df$Filters...exomes == "PASS" | GOI_gnomAD_df$Filters...genomes == "PASS",]
cat("\t\tFiltered variants:",length(GOI_gnomAD_df_filtered[,1]),"\n\n")

# remove synonymous variants
cat("\tFiltering out 'synonymous' variants")
GOI_gnomAD_df_filtered <- GOI_gnomAD_df_filtered[GOI_gnomAD_df_filtered$Annotation != "synonymous_variant",]
cat("\t\tFiltered variants:",length(GOI_gnomAD_df_filtered[,1]),"\n\n")

# remove intron variants
cat("\tFiltering out 'intron' variants (splice region kept)")
GOI_gnomAD_df_filtered <- GOI_gnomAD_df_filtered[GOI_gnomAD_df_filtered$Annotation != "intron_variant",]
cat("\t\tFiltered variants:",length(GOI_gnomAD_df_filtered[,1]),"\n\n")

# remove cases where Allele.Count = 0  ; not sure why these are returned at all ??
cat("\tFiltering out variants with Allele.Count of 0")
GOI_gnomAD_df_filtered <- GOI_gnomAD_df_filtered[GOI_gnomAD_df_filtered$Allele.Count != 0,]
cat("\t\tFiltered variants:",length(GOI_gnomAD_df_filtered[,1]),"\n\n")

save.image("troubleshooting_workspace.RData") #####################

cat("\tConverting coordinates from hg19 to hg38\n")
# this returns coordinates aligned to GRCh37/hg19; change to hg38
GOI_gnomAD_df_filtered$Position <- unlist(start(liftOver(GRanges(paste0('chr',GOI_CHR), IRanges(start = as.numeric(GOI_gnomAD_df_filtered$Position), width = 1)),Chain_19to38)))

save.image("troubleshooting_workspace.RData") #####################

#change ins and del format to match cBP (ExAC includes ref allele at start; e.g. ATT/A rather than TT/-)
cat("\tfixing ins and del format to match cBP variants\n\n")
InsDelIndex_gnomAD <- which((nchar(GOI_gnomAD_df_filtered$Reference)>1) | (nchar(GOI_gnomAD_df_filtered$Alternate)>1))
GOI_gnomAD_df_filtered[InsDelIndex_gnomAD,"Reference"] <- sapply(GOI_gnomAD_df_filtered$Reference[InsDelIndex_gnomAD], function(x) substring(x,2))
GOI_gnomAD_df_filtered[GOI_gnomAD_df_filtered$Reference == "","Reference"] <- "-"
GOI_gnomAD_df_filtered[InsDelIndex_gnomAD,"Alternate"] <- sapply(GOI_gnomAD_df_filtered$Alternate[InsDelIndex_gnomAD], function(x) substring(x,2))
GOI_gnomAD_df_filtered[GOI_gnomAD_df_filtered$Alternate == "","Alternate"] <- "-"

save.image("troubleshooting_workspace.RData") #####################


