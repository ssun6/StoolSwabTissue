## build command: docker build -f graphics.Dockerfile -t asorgen/graphics:v1 .

FROM r-base:3.6.0

#1.) set shell to bash
SHELL ["/bin/bash", "-c"]
ARG DEBIAN_FRONTEND=noninteractive

#2.) Install R Packages, ggplot2 and its error-prone dependencies
RUN Rscript -e "install.packages('testthat', dependencies=c('Depends', 'Imports') )" && \
    Rscript -e "install.packages('isoband', dependencies=c('Depends', 'Imports') )" && \
    Rscript -e "install.packages('ggplot2', dependencies=c('Depends', 'Imports') )" && \
	Rscript -e "install.packages('stringr', dependencies=c('Depends', 'Imports') )" 

#3.) Install more R Packages
RUN Rscript -e "install.packages('gridExtra', dependencies=c('Depends', 'Imports') )" && \
    Rscript -e "install.packages('ggsignif', dependencies=c('Depends', 'Imports') )" && \
    Rscript -e "install.packages('reshape2', dependencies=c('Depends', 'Imports') )" && \
    Rscript -e "install.packages('ggthemes', dependencies=c('Depends', 'Imports') )" && \
	Rscript -e "install.packages('ggrepel', dependencies=c('Depends', 'Imports') )"
RUN Rscript -e "install.packages('pheatmap', dependencies=c('Depends', 'Imports') )"
RUN Rscript -e "install.packages('nlme', dependencies=c('Depends', 'Imports') )"
RUN Rscript -e "install.packages('dplyr', dependencies=c('Depends', 'Imports') )"
RUN Rscript -e "install.packages('stringr', dependencies=c('Depends', 'Imports') )"

#4.) Install vegan R Packages
RUN Rscript -e "install.packages('vegan', dependencies=c('Depends', 'Imports') )"

#5.) Install more R Packages
#RUN Rscript -e "install.packages('tidyr', dependencies=c('Depends', 'Imports') )"
#RUN Rscript -e "install.packages('ggpubr', dependencies=c('Depends', 'Imports') )"
#RUN Rscript -e "install.packages('lme4', dependencies=c('Depends', 'Imports') )"
#RUN Rscript -e "install.packages('lmerTest', dependencies=c('Depends', 'Imports') )"
#RUN Rscript -e "install.packages('RGraphics', dependencies=c('Depends', 'Imports') )"


#6.) check that packages installed
RUN Rscript -e "library('ggplot2'); library('gridExtra'); library('reshape2'); library('ggrepel'); library('vegan'); "
RUN Rscript -e "library('ggthemes'); library('nlme'); library('pheatmap'); library('dplyr'); library('stringr'); "
#RUN Rscript -e "library('tidyr'); "
#RUN Rscript -e "library('ggpubr'); "
#RUN Rscript -e "library('lme4'); "
#RUN Rscript -e "library('lmerTest'); "
#RUN Rscript -e "library('RGraphics'); "


#7.) Cleanup
RUN	apt-get clean && \
	find / -name *python* | xargs rm -rf && \
	rm -rf /tmp/* && \
	rm -rf /usr/share/* && \
	rm -rf /var/cache/* && \
	rm -rf /var/lib/apt/lists/* && \
	rm -rf /var/log/*
