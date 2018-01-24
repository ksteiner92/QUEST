.SUFFIXES : .mod .o .F90

test  : $(QUEST_DIR)/$(DQMCLIB) test1.F90 test2.F90
	$(FC) $(FLAGS) -I$(QUEST_DIR)/SRC -o test1  test1.F90 $(QUEST_DIR)/$(DQMCLIB) $(LIB) 
	$(FC) $(FLAGS) -I$(QUEST_DIR)/SRC -o test2  test2.F90 $(QUEST_DIR)/$(DQMCLIB) $(LIB)

clean:
	rm -f test? 

