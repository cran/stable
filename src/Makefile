LIBNAME=stable
LD=ld
OBJS= stable.o	pstable.o

$(LIBNAME).so:	$(OBJS)
	@$(LD) $(SHLIBLDFLAGS) -o $(LIBNAME).so $(OBJS)

clean:
	@rm -f *.o *.so

realclean:
	@rm -f Makefile *.o *.so
