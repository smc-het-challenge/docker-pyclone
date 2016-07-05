FROM keyuan/docker-ccube
MAINTAINER "Ke Yuan" ke.yuan.09@gmail.com

COPY ./PyDP-0.2.2.tar.gz /home/pipeline/PyDP-0.2.2.tar.gz
COPY ./pyclone_ke.tar.gz /home/pipeline/pyclone_ke.tar.gz

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

COPY ./run_analysis_pyclone.R /home/pipeline/run_analysis_pyclone.R
RUN chmod +x /home/pipeline/run_analysis_pyclone.R 
