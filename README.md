# forest.din v0.1-beta.1.2

## Função para determinação da dinâmica de comunidades de espécies arbóreas                         

### Determinação das taxas de mortalidade, recrutamento, perda e ganho em área basal, mudanças líquidas e rotatividade, baseado em:

- KORNING, J.; BALSLEV, H. Growth and mortality of trees in Amazonian tropical rain forest in Ecuador. Journal of Vegetation Science,
v.5, n.1, p.77-86, 1994.
- OLIVEIRA FILHO, A. T. et a. Dinâmica da comunidade e populações arbóreas da borda e interior de um remanescente 
florestal na Serra da Mantiqueira, Minas Gerais, em um intervalo de cinco anos (1999-2004). 
Revista Brasileira de Botânica, v.30, n.1, p.149-161, 2007.
- SALAMI, B. et al. Influência de variáveis ambientais na dinâmica do componente arbóreo em um fragmento de Floresta
Ombrófila Mista em Lages, SC. Scientia Forestalis, v.42, n.102, p.197-207, 2014.
- SHEIL, D.; DAVID, BURSLEM, D. F. R. P.; ALDER, D. The interpretation and misinterpretation of mortality rate measures. Journal of Ecology, v.83, n.2, p.331-333, 1995.
- SHEIL, D.; JENNINGS, S.; SAVILL, P. Long-term permanent plot observations of vegetation dynamics in Budongo, a Ugandan rain forest. Journal of Tropical Ecology, v.16, n.6, p.865-882, 2000.

### Estimativa de Biomassa Acima do Solo baseado em Chave et al. (2014), considerando fórmula com as seguintes variáveis:                       

- D = Diâmetro (cm)
- WD = Densidade da madeira (g.cm-3)
- E = Estimativa de Stress Ambiental, baseado na coordenada geografica (coord)

* Referência:
- CHAVE et al. (2014) Improved allometric models to estimate the aboveground biomass of tropical trees, Global Change Biology, 20 (10), 3177-3190


#### Autor:  Pedro Higuchi                                   
 01/04/2017	
* Como citar

* HIGUCHI, P. forest.din: Função em linguagem de programação estatística R para determinação da dinâmica de comunidades de espécies arbóreas 2018. DOI: 10.5281/zenodo.1297702 Disponvel em https://github.com/higuchip/forest.din

* REJOU-MECHAIN, M.; TANGUY, A.; PIPONIOT, C.; CHAVE, J.; HERAULT, B. 	BIOMASS: Estimating Aboveground Biomass and Its Uncertainty in Tropical Forests. R package version 1.2. https://CRAN.R-project.org/package=BIOMASS	

													                           
# Observações:											                      
- a) O argumento x (planilha de dados) terá que conteras colunas Parcelas (identificação das parcelas),	Especie (id. espécies), DAP1 (DAP no ano 1) e  DAP2 (DAP no ano 2)   
- b) arquivo exemplo de entrada, disponível em https://raw.githubusercontent.com/higuchip/forest.din/master/dados_exemplo.csv
- c) O argumento t, representa o tempo entre inventários  
- d) Argumento coord deve ser do tipo c(long, lat), com valores em graus decimais
- e) Argumento add_wd representa um data.frame com valores de densidade da madeira (g.cm-3) formato com três colunas ("genus", "species", "wd"). Caso argumento add_wd não seja fornecido, a densidade da madeira será estimada com a função getWoodDensity do pacote BIOMASS, baseado em Zanne et al. Global wood density database. Dryad. Identifier: http://datadryad.org/handle/10255/dryad.235 (2009).
- f) Requer pacote BIOMASS


Caso tenha dúvidas, sugestões ou queira contribuir, entre em contato: higuchip@gmail.com

