if(rgplates:::getOS()=="linux") wd <- file.path(Sys.getenv("Dropbox"), "Software/biodome")
if(rgplates:::getOS()=="windows") wd <- file.path("D:/biodome")
if(rgplates:::getOS()=="osx") wd <- file.path("~/Desktop/biodome")
