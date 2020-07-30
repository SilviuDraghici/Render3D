CFLAGS = -O3 -g -lm -fopenmp
UOBJ = timer.o utils.o affineTransforms.o imageProcessor.o objects.o BVH.o meshes.o mappings.o ray.o camera.o buildscene.o rayTracer.o pathTracer.o Render3D.o
UTILS = utils/
RAY = rayTracer/
PTH = pathTracer/

SCENES = $(wildcard scenes/buildScene*.cpp)

Render3D: ${UOBJ}
	@echo Building Render3D !!!
	g++ ${CFLAGS} -o $@ ${UOBJ}
	@echo

OBJ = obj/

Render3D.o: Render3D.cpp  ${RAY}rayTracer.h
	@echo Building Render3D.o
	g++ ${CFLAGS} -c $^
	@echo

timer.o: ${UTILS}timer.cpp ${UTILS}timer.h
	@echo Building Timer
	g++ ${CFLAGS} -c $^
	@echo

utils.o: ${UTILS}utils.cpp ${UTILS}utils.h
	@echo Building Utils
	g++ ${CFLAGS} -c $^
	@echo

affineTransforms.o: ${UTILS}affineTransforms.cpp ${UTILS}affineTransforms.h
	@echo Building Affine Tranforms
	g++ ${CFLAGS} -c $^
	@echo

imageProcessor.o: ${UTILS}imageProcessor.cpp ${UTILS}imageProcessor.h
	@echo Building Image Processor
	g++ ${CFLAGS} -c $^
	@echo

objects.o: ${UTILS}objects.cpp ${UTILS}objects.h ${UTILS}color.h ${UTILS}mappings.h ${UTILS}ray.h ${UTILS}utils.h 
	@echo Building Objects
	g++ ${CFLAGS} -c $^
	@echo

BVH.o: ${UTILS}BVH.cpp ${UTILS}BVH.h  
	@echo Building BVH
	g++ ${CFLAGS} -c $^
	@echo

meshes.o: ${UTILS}meshes.cpp ${UTILS}meshes.h ${UTILS}color.h ${UTILS}mappings.h ${UTILS}ray.h ${UTILS}utils.h ${UTILS}objects.h ${UTILS}BVH.h
	@echo Building Meshes
	g++ ${CFLAGS} -c $^
	@echo

mappings.o: ${UTILS}mappings.cpp ${UTILS}mappings.h 
	@echo Building Mapping Funtions
	g++ ${CFLAGS} -c $^
	@echo

ray.o: ${UTILS}ray.cpp ${UTILS}ray.h ${UTILS}affineTransforms.h ${UTILS}objects.h ${UTILS}utils.h
	@echo Building Ray
	g++ ${CFLAGS} -c $^
	@echo

camera.o: ${UTILS}camera.cpp ${UTILS}camera.h
	@echo Building Camera
	g++ ${CFLAGS} -c $^
	@echo

buildscene.o: ${UTILS}buildscene.cpp ${UTILS}buildscene.h ${UTILS}ray.h ${UTILS}affineTransforms.h ${UTILS}objects.h ${UTILS}meshes.h ${UTILS}utils.h ${SCENES}
	@echo Building buildscene
	g++ ${CFLAGS} -c ${UTILS}buildscene.cpp ${UTILS}buildscene.h ${UTILS}ray.h ${UTILS}affineTransforms.h ${UTILS}objects.h ${UTILS}meshes.h ${UTILS}utils.h
	@echo

rayTracer.o: ${RAY}rayTracer.cpp ${RAY}rayTracer.h ${UTILS}buildscene.h ${UTILS}affineTransforms.h ${UTILS}imageProcessor.h ${UTILS}utils.h
	@echo Building Ray Tracer !
	g++ ${CFLAGS} -c $^
	@echo

pathTracer.o: ${PTH}pathTracer.cpp ${PTH}pathTracer.h ${UTILS}buildscene.h ${UTILS}affineTransforms.h ${UTILS}imageProcessor.h ${UTILS}utils.h
	@echo Building Path Tracer !
	g++ ${CFLAGS} -c $^
	@echo
	
clean:
	rm *.o
	rm utils/*.h.gch