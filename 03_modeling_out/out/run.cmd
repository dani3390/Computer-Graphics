.\bin\ymodel --scene tests\01_terrain\terrain.json --output outs\01_terrain\terrain.json --terrain object
.\bin\ymodel --scene tests\02_displacement\displacement.json --output outs\02_displacement\displacement.json --displacement object
.\bin\ymodel --scene tests\03_hair1\hair1.json --output outs\03_hair1\hair1.json --hairbase object --hair hair
.\bin\ymodel --scene tests\03_hair2\hair2.json --output outs\03_hair2\hair2.json --hairbase object --hair hair --hairlen 0.005 --hairstr 0
.\bin\ymodel --scene tests\03_hair3\hair3.json --output outs\03_hair3\hair3.json --hairbase object --hair hair --hairlen 0.005 --hairstr 0.01
.\bin\ymodel --scene tests\03_hair4\hair4.json --output outs\03_hair4\hair4.json --hairbase object --hair hair --hairlen 0.02 --hairstr 0.001 --hairgrav 0.0005 --hairstep 8
.\bin\ymodel --scene tests\04_grass\grass.json --output outs\04_grass\grass.json --grassbase object --grass grass


.\bin\yscene render outs\01_terrain\terrain.json --output out\01_terrain.jpg --samples 256 --resolution  720
.\bin\yscene render outs\02_displacement\displacement.json --output out\02_displacement.jpg --samples 256 --resolution  720
.\bin\yscene render outs\03_hair1\hair1.json --output out\03_hair1.jpg --samples 256 --resolution  720
.\bin\yscene render outs\03_hair2\hair2.json --output out\03_hair2.jpg --samples 256 --resolution  720
.\bin\yscene render outs\03_hair3\hair3.json --output out\03_hair3.jpg --samples 256 --resolution  720
.\bin\yscene render outs\03_hair4\hair4.json --output out\03_hair4.jpg --samples 256 --resolution  720
.\bin\yscene render outs\04_grass\grass.json --output out\04_grass.jpg --samples 256 --resolution  720 --bounces 128

:: Scene density
.\bin\ymodel --scene tests\03_hair1\hair1.json --output outs\03_hair1\hair1.json --hairbase object --hair hair --density 0.9
.\bin\ymodel --scene tests\03_hair2\hair2.json --output outs\03_hair2\hair2.json --hairbase object --hair hair --hairlen 0.005 --hairstr 0 --density 0.9


:: Render density
.\bin\yscene render outs\03_hair1\hair1.json --output out\density1.jpg --samples 256 --resolution  720
.\bin\yscene render outs\03_hair2\hair2.json --output out\density2.jpg --samples 256 --resolution  720

:: Extra2 u = 0
.\bin\ymodel --scene tests\02_displacement\displacement.json --output outs\02_displacement\displacement.json --displacement object --u 0 --v 0.05 --extra2

.\bin\yscene render outs\02_displacement\displacement.json --output out\extra2_v.jpg --samples 256 --resolution  720

:: Extra2 v = 0
.\bin\ymodel --scene tests\02_displacement\displacement.json --output outs\02_displacement\displacement.json --displacement object --u 0.5 --v 0 --extra2

.\bin\yscene render outs\02_displacement\displacement.json --output out\extra2_u.jpg --samples 256 --resolution  720

:: smoothvoronoi
.\bin\ymodel --scene tests\02_displacement\displacement.json --output outs\02_displacement\displacement.json --displacement object --u 1 --smooth

.\bin\yscene render outs\02_displacement\displacement.json --output out\smooth.jpg --samples 256 --resolution  720

:: voronoise
.\bin\ymodel --scene tests\02_displacement\displacement.json --output outs\02_displacement\displacement.json --displacement object --u 0 --v 1 --voronoise

.\bin\yscene render outs\02_displacement\displacement.json --output out\voronoise.jpg --samples 256 --resolution  720

:: cellnoise
.\bin\ymodel --scene tests\02_displacement\displacement.json --output outs\02_displacement\displacement.json --displacement object --u 0 --v 0 --cellnoise

.\bin\yscene render outs\02_displacement\displacement.json --output out\cellnoise.jpg --samples 256 --resolution  720




