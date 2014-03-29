#!/bin/bash

export VERSION=1.30
export GEPARDDIR=/home/jan/workspace/Gepard
export DISTDIR=/home/jan/workspace/Gepard/geparddist
export KEYSTORE=/home/jan/workspace/Gepard/keys/gepard

echo
echo "gepard.jar exported ?????"
echo
read

#      CREATE LIB DIR
mkdir $DISTDIR/lib

#      SIGN THE JAR
jarsigner -storepass plattformesdeka91 -keystore $KEYSTORE $DISTDIR/gepard.jar gepard 

#      BACKUP GEPARD.JAR
cp $DISTDIR/gepard.jar $DISTDIR/../ 

#      COPY ADDITIONAL JAR FILES
cp $GEPARDDIR/lib/jaxrpc.jar $DISTDIR/lib/
cp $GEPARDDIR/lib/saaj.jar $DISTDIR/lib/
cp $GEPARDDIR/lib/axis.jar $DISTDIR/lib/
cp $GEPARDDIR/lib/commons-logging-1.0.4.jar $DISTDIR/lib/
cp $GEPARDDIR/lib/commons-discovery-0.2.jar $DISTDIR/lib/
cp $GEPARDDIR/lib/wsdl4j-1.5.1.jar $DISTDIR/lib/
cp $DISTDIR/gepard.jar $DISTDIR/lib/

#      SIGN THE JARS
jarsigner -storepass plattformesdeka91 -keystore $KEYSTORE $DISTDIR/lib/gepard.jar gepard 
jarsigner -storepass plattformesdeka91 -keystore $KEYSTORE $DISTDIR/lib/jaxrpc.jar gepard 
jarsigner -storepass plattformesdeka91 -keystore $KEYSTORE $DISTDIR/lib/saaj.jar gepard 
jarsigner -storepass plattformesdeka91 -keystore $KEYSTORE $DISTDIR/lib/axis.jar gepard 
jarsigner -storepass plattformesdeka91 -keystore $KEYSTORE $DISTDIR/lib/commons-logging-1.0.4.jar gepard 
jarsigner -storepass plattformesdeka91 -keystore $KEYSTORE $DISTDIR/lib/commons-discovery-0.2.jar gepard 
jarsigner -storepass plattformesdeka91 -keystore $KEYSTORE $DISTDIR/lib/wsdl4j-1.5.1.jar gepard 

#      MOVE TO SUBDIR
cd $DISTDIR/

#      COPY ALL JAR TO THEIR OWN DIRECTORY
mkdir ../web
cp lib/*.jar ../web/

#     COPY ALL REMAINING FILES
cp $GEPARDDIR/stuff/gepard*.bat $DISTDIR
cp $GEPARDDIR/stuff/*.sh $DISTDIR
cp $GEPARDDIR/stuff/README.txt $DISTDIR
cp $GEPARDDIR/stuff/tutorial.html $DISTDIR
mkdir $DISTDIR/matrices
cp $GEPARDDIR/matrices/*.mat $DISTDIR/matrices/

#     MAKE .SH SCRIPTS EXECUTABLE 
chmod a+rx $DISTDIR/*.sh

#     DO PACKING
mkdir $DISTDIR/gepard-$VERSION
mv * gepard-$VERSION

#      PACK THE STUFF AS ZIP
zip -r -9 ../gepard-$VERSION.zip gepard-$VERSION 
#      PACK THE STUFF AS TAR.GZ
tar -pcz --file=../gepard-$VERSION.tar.gz gepard-$VERSION

# CLEAN UP
rm gepard-$VERSION gepard.jar -rf

mv ../gepard.jar .
mv ../gepard*.zip .
mv ../gepard*.tar.gz .
mv ../web .

#   COPY WEBSTART DESCRIPTOR FILES
cp $GEPARDDIR/stuff/*.jnlp $DISTDIR
