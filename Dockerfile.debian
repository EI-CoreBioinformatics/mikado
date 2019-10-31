FROM conda/miniconda3

LABEL description="Docker file for the Mikado pipeline"
WORKDIR /
RUN apt update && apt-get install -qq -y git wget zlib1g-dev g++ && apt clean
COPY environment.yml /
RUN conda env update --prune -n base -f /environment.yml && conda clean -afy
RUN git clone https://github.com/EI-CoreBioinformatics/mikado.git /usr/local/src/mikado
WORKDIR /usr/local/src/mikado
RUN python setup.py bdist_wheel && pip install dist/*whl
RUN chash=$(git log | head -n1 | cut -f 2 -d " ") && echo -e "#!/bin/bash\necho ${chash}" > /usr/local/bin/show_commit_hash && chmod 775 /usr/local/bin/show_commit_hash
RUN chmod -R 775 /usr/local/src/mikado/util/*
RUN for TOOL in /usr/local/src/mikado/util/*; do script=$(basename ${TOOL}) && cp ${TOOL} /usr/local/bin/${script}; done
WORKDIR /usr/local/src/
RUN rm -rf mikado
CMD mikado
CMD daijin