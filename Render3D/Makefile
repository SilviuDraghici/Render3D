UOBJ = Render3D.o rayTracer.o utils.o affineTransforms.o imageProcessor.o objects.o mappings.o ray.o  
UTILS = utils/
RAY = rayTracer/

Render3D: ${UOBJ}
	@echo Building Render3D !!!
	g++ -O4 -g -lm -o $@ ${UOBJ}

OBJ = obj/

Render3D.o: Render3D.c  ${RAY}rayTracer.h
	@echo Building Render3D.o
	g++ -O4 -g -lm -c $^

utils.o: ${UTILS}utils.c ${UTILS}utils.h
	@echo Building Utils
	g++ -O4 -g -lm -c $^

affineTransforms.o: ${UTILS}affineTransforms.c ${UTILS}affineTransforms.h ${UTILS}svdDynamic.c ${UTILS}svdDynamic.h 
	@echo Building Affine Tranforms
	g++ -O4 -g -lm -c $^

imageProcessor.o: ${UTILS}imageProcessor.c ${UTILS}imageProcessor.h
	@echo Building Image Processor
	g++ -O4 -g -lm -c $^

objects.o: ${UTILS}objects.c ${UTILS}objects.h ${UTILS}color.h ${UTILS}mappings.h ${UTILS}ray.h ${UTILS}utils.h 
	@echo Building Objects
	g++ -O4 -g -lm -c $^

mappings.o: ${UTILS}mappings.c ${UTILS}mappings.h 
	@echo Building Mapping Funtions
	g++ -O4 -g -lm -c $^

ray.o: ${UTILS}ray.c ${UTILS}ray.h ${UTILS}affineTransforms.h ${UTILS}objects.h ${UTILS}utils.h
	@echo Building Ray
	g++ -O4 -g -lm -c $^

rayTracer.o: ${RAY}rayTracer.c ${RAY}rayTracer.h ${UTILS}affineTransforms.h ${UTILS}imageProcessor.h
	@echo Building Ray Tracer !
	g++ -O4 -g -lm -c $^

clean:
	rm *.o
	rm utils/*.h.gch