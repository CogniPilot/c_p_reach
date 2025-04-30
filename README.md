C-P Reach
---------------

Cyber Physical System Reachability Analyzer.


# Running as a docker container 


Then you can run it with a make target as well

```
make run 
```

Then connect to the webpage at `localhost:8000` 


## Building 
Install docker: 

Install build-essential:
```
apt install build-essential
```

Use the provided make target to build the container:

```
make build
```
## Running from docker 



To run from docker you can use a make target

```
make run 
```

Then connect to the webpage at `localhost:8000` 


# Running locally on your machine 
## Installation

* Install [Poetry](https://python-poetry.org/docs/)

* Install CDD Lib
```bash
sudo apt install libcdd-dev
```

* Build Poetry Package

```bash
git clone git@github.com:cognipilot/c_p_reach
cd c_p_reach
poetry install
```

## Running Commands

First start a poetry shell in the poetry virtual environment.

```bash
poetry shell
```
Next, call the reachability tool.
```bash
c_p_reach
```

# Running on a k8s cluster

# Updating the runtime parameters 

Update `chart/values.yml` to point to your docker registry and k8s cluster

```
# Cluster host is the dns address of your k8s ingress
cluster_host: whatever.com 
# Exposed url is the subnet at which to expose the tool the url will be https://<exposed_url>.<cluster_host>
# So for this example the tool will be reachable at https://cp-reach.whatever.com
exposed_url: cp-reach    


images:
  # This is the url for your docker registry where the cpreach image is located
  url: myDockerRegistry.com:5050
  pullPolicy: Always
  # This is the login credentials for your docker registry 
  username: username
  password: password
  credentialName: registry-creds

gui:
  # This is the full image path in the docker registry to the cp-reach image
  # so in this example the cp-reach image is at myDockerRegistry.com:5050/tools/cp-reach/gui:latest
  image: tools/cp-reach/gui:latest
  # dont change port
  port: 8501


```


Generate a `.kubecfg` file for your k8s cluster and save it at `./kubectx` (do not check into git!)



# Roadmap

## Working
* Multirotor 3D support
* Support Rover

## TODO
* Auto-read gains from CasADi model
* Take vehicle trajectory as input
