Plasma Physics simulation

by
Elias Greil, Sebastian Riepl

How to use CWS:

    Linux prerequisite [preferred]:
    docker, docker-compose
    
    Windows:
    DockerHub Desktop


- either build the container or pull it from docker hub [recommended] 
    docker pull sebobo233/plasma:1.0

- rename the image to an apropriate name, e.g. CWS
    docker tag sebobo233/plasma:1.0 CWS:1.0

- start the container with:
    on Linux use the docker-compose file
        - go to the file location where docker-compose is saved
        - type "docker compose up -d"
        - look for the login key of the jupyternotebook in the log by using: "docker logs CWS"

    on Windows use docker desktop for downloading and starting it

- connect to jupyter server by typing in browser on the same machine:
    127.0.0.1:8080
    login with key from log, see step above   



If you need help contact me on Discord:
We are in the coil optimization channel [PN to user: sebobo]




important sources:

https://quasr.flatironinstitute.org/model/0021326

https://github.com/itpplasma/ALPES/tree/main

https://github.com/hiddenSymmetries/simsopt

https://simsopt.readthedocs.io/en/latest/overview.html


