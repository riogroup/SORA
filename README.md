# occultation
Nome Sugerido / Suggested name (Bruno Morgado)
# Stellar Occultation Reduction and Analysis (SORA)

SORA significa "céu" em japonês / SORA means "sky" in Japanese

# Descrição do programa (English follows)

Programa que atualiza e automatiza os programas do Bruno Sicardy utilizados para redução de dados de observações de ocultaçes estelares.

Aqui foram criadas e adaptadas as rotinas baseadas (e não 'traduzidas') nos programas* e subrotinas:
  - positionv
  - ephem_planete
  - fit_d2_ksi_eta
  - diam
  - bar
  - polyfit

*Todos os programas em fortran foram criados pelo Bruno Sicardy. As equações e a programação está (algumas vezes) documentada nos próprios pogramas, principalmente as baseadas em rotinas do "Numerical Recipes", e será reproduzida aqui conforme desenvolvimento do programa.

O novo programa é dividido em subrotinas, pacotes e classes de modo a facilitar modificações e atualizações em módulos específicos, sem precisar realizar grandes alterações no programa principal.

O programa também é pensado em adaptações para uma integração com o portal TNO do LIneA.

===
Incluir uma descrição do que o programa faz
Incluir o que o programa usa como INPUT e OUTPUT

**********

# Program description

Program to update and automate Bruno Sicardy's routines used for data reduction obtained in stellar occultations observations.

Here we create and adapted routines based (and not 'translated') on fortran programs* and subroutines:
  - positionv
  - ephem_planete
  - fit_d2_ksi_eta
  - diam
  - bar
  - polyfit
  
  *All fortran programs were created by Bruno Sicardy. Equations and programming is documented (some times) inside each code, mainly the ones based on routines from "Numerical Recipes", and it will be reproduced here as we develop the prgram.
  
The new program is divided in subroutines, packages and classes, in order to facilitade future modifications and new updates in specific modules, without the need to perform big changes in the main code.
  
This program was also developed in order to make easy adaptations for an integration with TNO portan at LIneA.


===
Include a description about what does the program do
Incluir what does it uses as INPUT and OUTPUT
