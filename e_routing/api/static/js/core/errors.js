/**
 * Shared application error presentation.
 */

function showError(message) {
    $('alert').hidden = false;
    $('alert').textContent = message
}

function clearError() {
    $('alert').hidden = true;
    $('alert').textContent = ''
}