CC=g++
SRC=tfserving_utils.cpp tfs_messages.cpp
OBJ=tfserving_utils.o tfs_messages.o
LIBS= -lcurl
default: all
all: $(OBJ)
	$(CC) -o exe $(OBJ) $(LIBS)
clean:
	rm  *.o
	rm z
tfserving_utils.o: tfserving_utils.cpp
	$(CC) $(CFLAGS) -c tfserving_utils.cpp -o tfserving_utils.o $(LIBS) 
tfs_messages.o: tfs_messages.cpp
	$(CC) $(CFLAGS) -c tfs_messages.cpp -o tfs_messages.o $(LIBS) 
