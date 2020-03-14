---
category: Core
path: '/Core/:id'
title: 'core.sb_programs'

layout: nil
---
The **core.sb_programs** contains methods for parsing command strings and using the appropriate container environment. The interaction between workflows and containers is carried out through these methods.  
<br>
## *class* core.sb_programs**Run**(command, path, docker_image)  
#### The Run object performs the function of preparing the command string to be run in a container and selecting the appropriate container environment.  
<br>
#### The Run object has several initialization parameters needed to function.
#### The *command* parameter is a string that provides the command that will be executed in the container.
#### The *path* parameter is a python dictionary where the key the path on the host file system and the value is the associated mounting location in the container environment eg: *{"/path/outside":"/path/incontainer"}*.
#### The *docker_image* parameter provides the name of the container environment that should be used according to the *docker_config.json* file. Note, that when using singularity the images from Docker are used to create the Singularity environment so the *docker_image* parameter still applies.
<br>
### **run()**  
#### Once the Run object has been created the command can be executed by calling the *run()* function on the object. This calls the appropriate container engine and begins executing the command stored in the object.
<pre><code class="python">from core.sb_programs import Run

command = "echo Hello World"
path = {"/home":"/container/home"}
container = "spades"

spades = Run(command=command,path=path,docker_image=container)
spades.run()
</code></pre>
<br>
## *class* core.sb_programs**Run_multi**(command_list, path, docker_image)  
#### The Run_multi object performs a similar function as Run by preparing a list of command strings to be run in multiple containers simultaneously and selecting the appropriate container environment.  
<br>
#### The Run object has several initialization parameters needed to function.
#### The *command_list* parameter is a list of strings that provides the commands that will each be executed in a separate but identical container.
#### The *path* parameter is a python dictionary where the key the path on the host file system and the value is the associated mounting location in the container environments eg: *{"/path/outside":"/path/incontainer"}*.
#### The *docker_image* parameter provides the name of the container environment that should be used according to the *docker_config.json* file. Note, that when using singularity the images from Docker are used to create the Singularity environment so the *docker_image* parameter still applies.
<br>
### **run(jobs)**  
#### Once the Run object has been created the commands can be executed by calling the *run()* function on the object. The *jobs* parameter specifies how many concurrent containers to run at once.
#### **Caution:** this is independent of how many CPUs are used by the program. Containerization by default makes all host resources available to the container. If the default nature of the program is to use all available CPUs then using the Run_multi object could result in attempting use more resources than available. This can result in a large reduction in system performance. If using the multiprocessing approach be sure to control the amount of system resources that are being used by the application in concert with the number of concurrent jobs.
<pre><code class="python">from core.sb_programs import Run_multi

commands = ["echo Hello World One","echo Hello World Two","echo Hello World three"]
path = {"/home":"/container/home"}
container = "spades"

spades = Run_multi(command_list=commands,path=path,docker_image=container)
spades.run(2)
</code></pre>
