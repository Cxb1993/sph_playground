docker-compose -f ./develop/devcompose.yml up -d
docker exec -it \
	-e LINES=$LINES \
	-e COLUMNS=$COLUMNS \
	-e TERM=$TERM \
	develop_fortran_1 \
	zsh
