if(rgplates:::getOS()=="linux") wd <- file.path(Sys.getenv("Dropbox"), "Software/orange")
if(rgplates:::getOS()=="windows") wd <- file.path("D:/orange")
if(rgplates:::getOS()=="osx") wd <- file.path("~/Desktop/orange")
