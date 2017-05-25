FROM  debian:9

RUN   cd ~ && \
      ulimit -n 1000 && \
      apt-get update && \
      apt-get install -y valgrind make gfortran zsh git-core wget && \
      wget https://github.com/robbyrussell/oh-my-zsh/raw/master/tools/install.sh -O - | zsh
