.SHELL=/bin/bash

.PHONY=clean
HTSJDK=${HOME}/src/htsjdk/build/libs/htsjdk-4.2.0-6-g6eb0607-SNAPSHOT.jar
TESTNG=${HOME}/src/jvarkit/lib/com/beust/jcommander/1.82/jcommander-1.82.jar:${HOME}/src/jvarkit/lib/org/testng/testng/6.14.3/testng-6.14.3.jar

ifeq (${JAVA_HOME},)
$(error $${JAVA_HOME} is not defined)
endif

all: bcf.jar testng.xml bcf.jar testdata/test01.bcf.csi
	java -cp ${TESTNG}:${HTSJDK}:bcf.jar org.testng.TestNG testng.xml

CLASSES=BCFCodec BCFFileReader BCFIterator

bcf.jar : $(addsuffix .java,$(addprefix ./src/main/java/com/github/lindenb/jvarkit/variant/bcf/,$(CLASSES) BCFTypedData)) \
	$(addsuffix Test.java,$(addprefix ./src/test/java/com/github/lindenb/jvarkit/variant/bcf/,$(CLASSES)))
	rm -rf tmp
	mkdir -p tmp
	javac -cp $(HTSJDK):$(TESTNG) -d tmp $^
	jar cvf $@ -C tmp .
	rm -rf tmp

testdata/test01.bcf.csi : testdata/test01.bcf
	bcftools index -f $<

testdata/test01.bcf : testdata/test01.vcf
	bcftools view -O b -o $@ $<
