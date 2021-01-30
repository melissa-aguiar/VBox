def main():
    f0 = open("textos.txt", "w+")
    f0.write("Ma = [Ma;")
    for i in range(64):
        f0.write(" Ma" + str(i+1) + "(i,:); ")
    f0.write("];")
    f0.close()

if __name__ == "__main__":
    main()