## build command: docker build -f tidy.Dockerfile -t asorgen/tidy:v1 .

FROM rocker/tidyverse:3.6.3

#1.) set shell to bash
SHELL ["/bin/bash", "-c"]
ARG DEBIAN_FRONTEND=noninteractive


#5.) Install more R Packages
RUN Rscript -e "install.packages('tidyr', dependencies=c('Depends', 'Imports') )"
RUN Rscript -e "install.packages('ggpubr', dependencies=c('Depends', 'Imports') )"



#6.) check that packages installed
RUN Rscript -e "library('ggpubr'); "
RUN Rscript -e "library('tidyr'); "

#7.) Cleanup
RUN	apt-get clean && \
	find / -name *python* | xargs rm -rf && \
	rm -rf /tmp/* && \
	rm -rf /usr/share/* && \
	rm -rf /var/cache/* && \
	rm -rf /var/lib/apt/lists/* && \
	rm -rf /var/log/*