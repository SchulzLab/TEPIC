ggplot<-require("gplots")
ggplot2<-require("ggplot2")
glmnet<-require("glmnet")
doMC<-require("doMC")
reshape2<-require("reshape2")
gridExtra<-require("gridExtra")

if (ggplot){
	print("gplots is already installed")
}else{
	print("Installing gplots")
	install.packages("gplots",repos="https://ftp.gwdg.de/pub/misc/cran/")
}

if (ggplot2){
	print("ggplot2 is already installed")
}else{
	print("Installing ggplot2")
	install.packages("ggplot2",repos="https://ftp.gwdg.de/pub/misc/cran/")
}

if (glmnet){
	print("glmnet is already installed")
}else{
	print("Installing glmnet")
	install.packages("glmnet",repos="https://ftp.gwdg.de/pub/misc/cran/")
}

if (doMC){
	print("doMC is already installed")
}else{
	print("Installing doMC")
	install.packages("doMC",repos="https://ftp.gwdg.de/pub/misc/cran/")
}

if (reshape2){
	print("reshape2 is already installed")
}else{
	print("Installing reshape2")
	install.packages("reshape2",repos="https://ftp.gwdg.de/pub/misc/cran/")
}

if (gridExtra){
	print("gridExtra is already installed")
}else{
	print("Installing gridExtra")
	install.packages("gridExtra",repos="https://ftp.gwdg.de/pub/misc/cran/")
}

