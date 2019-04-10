FROM conda/miniconda3

LABEL description="Docker file for the Mikado pipeline"
WORKDIR /
RUN apt update && apt-get install -qq -y git wget zlib1g-dev g++ && apt clean
COPY environment.yml /
RUN conda env update -n base -f /environment.yml
RUN conda env list
RUN git clone https://github.com/EI-CoreBioinformatics/mikado.git /usr/local/src/mikado
WORKDIR /usr/local/src/mikado
RUN python setup.py bdist_wheel && pip install dist/*whl
RUN echo -e "#!/bin/bash\ncd /usr/local/src/mikado;\ngit log | head -n1 | cut -f 2 -d ' '" > /usr/local/bin/show_commit_hash && chmod 775 /usr/local/bin/show_commit_hash
CMD mikado
CMD daijin