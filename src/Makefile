HTML=lcov -t test -o rep.info -c -d ./
CFLAGS=--std=c++17 -lstdc++ -lm
TESTFLAGS=-lgtest -lgcov
GCOVFLAGS=--coverage
OS = $(shell uname)

all: clean gcov_report

s21_matrix_oop.a: clean
	gcc $(GCOVFLAGS) -c s21_matrix.cc
	ar rc s21_matrix_oop.a s21_matrix.o
	ranlib s21_matrix_oop.a

test: s21_matrix_oop.a
	gcc test.cc s21_matrix_oop.a $(CFLAGS) $(TESTFLAGS) -o test
	./test

gcov_report: test
	$(HTML)
	genhtml -o report rep.info
	open report/index.html

clang:
	cp ../materials/linters/.clang-format .clang-format
	clang-format -i *.cc *.h
	rm -rf .clang-format

clang_review:
	cp ../materials/linters/.clang-format .clang-format
	clang-format -n *.cc *.h
	rm -rf .clang-format

clean:
	rm -rf *.o *.a *.so *.gcda *.gcno *.gch rep.info *.html *.css test report *.txt *.dSYM

check: test
	cppcheck --enable=all --suppress=missingIncludeSystem --inconclusive --check-config s21_matrix.cc *.h
ifeq ($(OS), Darwin)
	leaks --atExit -- test
else
	CK_FORK=no valgrind --vgdb=no --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --log-file=RESULT_VALGRIND.txt ./test
endif