# LightSense-Calib: Calibração Espectral do Sensor AS7341

## 📌 Visão Geral do Projeto
Este projeto faz parte do meu **Trabalho de Conclusão de Curso (TCC) na UFSM**, focado na calibração de um sensor de luz multicanal, especificamente o **OSRAM AS7341**. Mais detalhes podem ser encontrados no meu TCC: [Repositório UFSM](https://repositorio.ufsm.br/handle/1/32820).

## 🎯 Objetivo
O objetivo é calibrar o sensor AS7341 comparando os seus dados espectrais reconstruídos com espectros de referência de alta fidelidade, obtidos por meio de uma **esfera integradora**. Os algoritmos de calibração visam encontrar o melhor conjunto de coeficientes que minimizem o erro entre o espectro reconstruído pelo sensor e o espectro de referência.

## 🔬 Processo de Calibração
O projeto possui três scripts principais no MATLAB:

### **1️⃣ channels_calib_no_diff.m**
- Este script calibra cada um dos **9 canais espectrais** do sensor AS7341.
- Gera **coeficientes de calibração** para uso posterior.
- Produz **gráficos e imagens** detalhando cada etapa do processo.

  ![Calibração sem difusor](saved_images/CalibrationIdealCoefficients_RB_100.png)

### **2️⃣ channels_calib_with_diff.m**
- Utiliza os coeficientes da calibração anterior.
- Faz a calibração do sensor **com difusor**, essencial para garantir medições precisas de luz.
- Gera os **coeficientes finais** considerando o difusor.

  ![Comparação métodos de calibração](saved_images/Rsquared_eachMethod_eachTrial_diffusers.png)

### **3️⃣ sensor_processor.m**
- Recebe os dados do sensor e os processa **com ou sem calibração**.
- Exibe medidas como **Iluminância e PPFD (Densidade de Fluxo de Fótons Fotossintéticos)**.

  ![Espectro reconstruído RGB](saved_images/SPD_reconstructed_methodsRGB_100.pdf)

## 🏆 Metodologia de Calibração
Os dados de calibração são obtidos por meio da comparação entre:
1. O **espectro reconstruído pelo sensor** AS7341.
2. O **espectro de referência** da mesma fonte de luz, obtido dentro de uma **esfera integradora** (equipamento de alta precisão).

Os algoritmos otimizam os coeficientes de calibração para minimizar os erros entre essas duas medições.

### 📊 Fluxo de calibração do sensor AS7341  
A imagem abaixo ilustra o processo de calibração do sensor, destacando as etapas de otimização dos coeficientes, validação via \(R^2\) e comparação dos resultados com uma esfera integradora:

## 📜 Mais Informações
Para uma explicação detalhada da metodologia, resultados e implementação, consulte o meu TCC: [Repositório UFSM](https://repositorio.ufsm.br/handle/1/32820).
