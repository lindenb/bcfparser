.SHELL=/bin/bash

.PHONY=clean test all jar javadoc

CLASSES=BCFCodec BCFFileReader BCFIterator

ifeq (${JAVA_HOME},)
$(error $${JAVA_HOME} is not defined)
endif



all: jar javadoc test

test : bcf.jar testng.xml bcf.jar lib/htsjdk.jar lib/jcommander.jar lib/testng.jar lib/jexl.jar
	java -cp lib/jcommander.jar:lib/testng.jar:lib/htsjdk.jar:bcf.jar:lib/jexl.jar org.testng.TestNG testng.xml



jar : bcf.jar

bcf.jar : $(addsuffix .java,$(addprefix ./src/main/java/com/github/lindenb/jvarkit/variant/bcf/,$(CLASSES) BCFTypedData)) \
	$(addsuffix Test.java,$(addprefix ./src/test/java/com/github/lindenb/jvarkit/variant/bcf/,$(CLASSES))) \
	lib/htsjdk.jar lib/jcommander.jar lib/testng.jar
	rm -rf tmp
	mkdir -p tmp
	javac -cp lib/htsjdk.jar:lib/jcommander.jar:lib/testng.jar -d tmp $(filter %.java,$^)
	jar cvf $@ -C tmp .
	rm -rf tmp

javadoc : $(addsuffix .java,$(addprefix ./src/main/java/com/github/lindenb/jvarkit/variant/bcf/,$(CLASSES))) lib/htsjdk.jar
	mkdir -p doc
	javadoc -cp lib/htsjdk.jar -d doc -sourcepath src/main/java com.github.lindenb.jvarkit.variant.bcf

lib/htsjdk.jar:
	mkdir -p $(dir $@)
	wget -O $@ "https://repo1.maven.org/maven2/com/github/samtools/htsjdk/4.3.0/htsjdk-4.3.0.jar"
	
lib/jcommander.jar:
	mkdir -p $(dir $@)
	wget -O $@ "https://repo1.maven.org/maven2/org/jcommander/jcommander/2.0/jcommander-2.0.jar"

lib/testng.jar:
	mkdir -p $(dir $@)
	wget -O $@ "https://repo1.maven.org/maven2/org/testng/testng/6.14.3/testng-6.14.3.jar"		

lib/jexl.jar:
	mkdir -p $(dir $@)
	wget -O $@ "https://repo1.maven.org/maven2/org/apache/commons/commons-jexl/2.1.1/commons-jexl-2.1.1.jar"

clean:
	rm -rf lib bcf.jar test-output

