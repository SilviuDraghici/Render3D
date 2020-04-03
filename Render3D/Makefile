UOBJ = utils.o affineTransforms.o imageProcessor.o objects.o mappings.o ray.o camera.o buildscene.o rayTracer.o pathTracer.o Render3D.o
UTILS = utils/
RAY = rayTracer/
PTH = pathTracer/

Render3D: ${UOBJ}
	@echo Building Render3D !!!
	g++ -O4 -g -lm -o $@ ${UOBJ}
	@echo

OBJ = obj/

Render3D.o: Render3D.c  ${RAY}rayTracer.h
	@echo Building Render3D.o
	g++ -O4 -g -lm -c $^
	@echo

utils.o: ${UTILS}utils.c ${UTILS}utils.h
	@echo Building Utils
	g++ -O4 -g -lm -c $^
	@echo

affineTransforms.o: ${UTILS}affineTransforms.c ${UTILS}affineTransforms.h ${UTILS}svdDynamic.c ${UTILS}svdDynamic.h 
	@echo Building Affine Tranforms
	g++ -O4 -g -lm -c $^
	@echo

imageProcessor.o: ${UTILS}imageProcessor.c ${UTILS}imageProcessor.h
	@echo Building Image Processor
	g++ -O4 -g -lm -c $^
	@echo

objects.o: ${UTILS}objects.c ${UTILS}objects.h ${UTILS}color.h ${UTILS}mappings.h ${UTILS}ray.h ${UTILS}utils.h 
	@echo Building Objects
	g++ -O4 -g -lm -c $^
	@echo

mappings.o: ${UTILS}mappings.c ${UTILS}mappings.h 
	@echo Building Mapping Funtions
	g++ -O4 -g -lm -c $^
	@echo

ray.o: ${UTILS}ray.c ${UTILS}ray.h ${UTILS}affineTransforms.h ${UTILS}objects.h ${UTILS}utils.h
	@echo Building Ray
	g++ -O4 -g -lm -c $^
	@echo

camera.o: ${UTILS}camera.c ${UTILS}camera.h
	@echo Building Camera
	g++ -O4 -g -lm -c $^
	@echo

buildscene.o: ${UTILS}buildscene.c ${UTILS}buildscene.h ${UTILS}ray.h ${UTILS}affineTransforms.h ${UTILS}objects.h ${UTILS}utils.h   
	@echo Building buildscene
	g++ -O4 -g -lm -c $^
	@echo

rayTracer.o: ${RAY}rayTracer.c ${RAY}rayTracer.h ${UTILS}buildscene.h ${UTILS}affineTransforms.h ${UTILS}imageProcessor.h ${UTILS}utils.h
	@echo Building Ray Tracer !
	g++ -O4 -g -lm -c $^
	@echo

pathTracer.o: ${PTH}pathTracer.c ${PTH}pathTracer.h ${UTILS}buildscene.h ${UTILS}affineTransforms.h ${UTILS}imageProcessor.h ${UTILS}utils.h
	@echo Building Path Tracer !
	g++ -O4 -g -lm -c $^
	@echo
	
clean:
	rm *.o
	rm utils/*.h.gch