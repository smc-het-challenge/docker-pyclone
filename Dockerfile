FROM keyuan/docker-ccube
MAINTAINER "Ke Yuan" ke.yuan.09@gmail.com

COPY ./PyDP-0.2.2.tar.gz /home/pipeline/PyDP-0.2.2.tar.gz
COPY ./pyclone_ke.tar.gz /home/pipeline/pyclone_ke.tar.gz

RUN apt-get update && apt-get upgrade -y python-pip

RUN pip install pyyaml \
    && cd /home/pipeline/ \
    && tar xvfz /home/pipeline/PyDP-0.2.2.tar.gz \
    && cd /home/pipeline/PyDP-0.2.2/ \
    && python setup.py install \
    && cd /home/pipeline/ \
    && tar xvfz /home/pipeline/pyclone_ke.tar.gz \
    && cd /home/pipeline/pyclone/ \
    && python setup.py install \
    && cd /home/

COPY ./create_ccfclust_inputs.py /home/pipeline/create_ccfclust_inputs.py
COPY ./run_analysis_pyclone.R /home/pipeline/run_analysis_pyclone.R

RUN chmod +x /home/pipeline/create_ccfclust_inputs.py \
    && chmod +x /home/pipeline/run_analysis_pyclone.R \
    && chmod +w /usr/local/lib/R/site-library
    
ENV PATH=/home/pipeline:$PATH
