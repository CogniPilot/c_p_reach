
.PHONY: build run deploy
build: 
	docker build -t cp-reach . 
run:
	docker run -it -p 8501:8501 cp-reach 
deploy:
	helm install cp-reach ./chart --kubeconfig ./kubectx --namespace cp-reach --create-namespace