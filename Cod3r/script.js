document.addEventListener('DOMContentLoaded', (event) => {
    const display = document.getElementById('display');
    const buttons = Array.from(document.getElementsByTagName('button'));
    
    buttons.map(button => {
        button.addEventListener('click', (e) => {
            const buttonText = e.target.innerText;
            
            if (buttonText === 'C') {
                display.value = '';
            } else if (buttonText === '=') {
                try {
                    display.value = eval(display.value);
                } catch (error) {
                    display.value = 'Erro';
                }
            } else {
                display.value += buttonText;
            }
        });
    });
});
