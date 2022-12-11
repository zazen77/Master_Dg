with open("Tunnel_excavation.txt", "w") as file:
    # Case of modeling
    # case 1 : flat Terrain and no buildings 
    # case 2 : flat Terrain and buildings
    # case 3 : real Terrain and no buildigns
    # case 4 : real Terrain and buildings
    
    case = 1

    if case == 1:
        case = 'X'
    elif case == 2: 
        case = 'B'
    elif case == 3: 
        case = 'G'
    else:
        case = 'BG'

    step = 0
    for i in range(1, 28): 
        step = step + 1 
        Mstep = 28 - step
        id = step + 10

        file.write("; Step" + str(step) + "\n")
        file.write("zone cmodel assign null range group 'ZG_M"+ str(Mstep) +"'\n")
        file.write("zone cmodel assign null range group 'ZG_N"+ str(step) +"'\n")
        file.write("struct shell create by-zone-face internal id " + str(id) + " range group 'EF_ML" + str(Mstep) + "' or 'EF_NL" + str(step) + "'\n")
        file.write("struct shell property isotropic (10.5e9,0.25) ...\n")
        file.write("                        thickness 0.3 ...\n")
        file.write("                        density 2500 ...\n")
        file.write("                        range id" + str(id) + "\n")
        file.write("model solve convergence 10\n")
        file.write("; Save\n")
        file.write("model save 'ZG_"+ str(step) +"'\n")
        file.write("plot 'disp' current\n")
        file.write("plot export bitmap dpi 300 size 1920 1080 filename '"+ str(step) +"_disp_plot.png'\n")
        file.write(";Mail send\n")
        file.write("program mail add to 'jaeeuncho.mc@gmail.com'\n")
        file.write("program mail subject 'step "+ str(step) +"'\n")
        file.write("program mail add attachment '"+ str(step) +"_disp_plot.png'\n")
        file.write(";----------------------------------------------------------------\n")
