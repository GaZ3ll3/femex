include ../../Makefile.in

coinipopt = Ipopt-3.11.9

url = http://www.coin-or.org/download/source/Ipopt/$(coinipopt).tgz

target = ipopt.tgz


all:	
	@echo " "
	@echo " "
	@echo "****** Building ipopt for femex ******"
	@echo " "
	
	@echo "0. Clean up previous downloaded tar..."
	@# rm -f $(target)
	@echo " "
	
	@echo "1. Downloading Ipopt from official website..."
	@echo " "

	@echo '        wget $(url)'
	@echo " "
	@wget $(url)
	@mv $(coinipopt).tgz $(target)

	
	@echo "2. Uncompressing the tarball..."
	@gunzip -f ipopt.tgz
	@echo " "
	@tar xf ipopt.tar
	@rm ipopt.tar
	
	
	@echo " "
	@echo "3. Downloading Third Party packages..."
	
	@cd $(coinipopt)/ThirdParty/Blas && ./get.Blas
	@cd $(coinipopt)/ThirdParty/Lapack && ./get.Lapack
	@cd $(coinipopt)/ThirdParty/Metis && ./get.Metis
	@cd $(coinipopt)/ThirdParty/Mumps && ./get.Mumps
	@cd $(coinipopt)/ThirdParty/ASL && ./get.ASL
	
	@echo " "
	@echo "Done downloading all Third Party packages."
	@echo " "
	
	@echo " "
	@echo "4. Compiling Ipopt..."
	@cd $(coinipopt) && ./configure CC=gcc F77=gfortran-4.4 CXX=g++ --disable-shared --with-pic\
	 --with-blas=BUILD --with-lapack=BUILD --with-matlab-home=$(MATLAB_ROOT) coin_skip_warn_cxxflags=yes
	@cd $(coinipopt) && make -j4
	@cd $(coinipopt) && make install
 	
	@echo " "
	@echo "5. Compiling with MATLAB_ROOT as following:"
	@echo " "
	@echo '        $(MATLAB_ROOT)'
	
	@rm -rf ipopt
	@mkdir ipopt
	@cd $(coinipopt)/Ipopt/contrib/MatlabInterface/src && make && mv ipopt.mexa64 ../../../../../ipopt/
	
	rm -rf $(coinipopt)
	
	@echo "Done......"
	
	
	

 	
 	
	

	
	
