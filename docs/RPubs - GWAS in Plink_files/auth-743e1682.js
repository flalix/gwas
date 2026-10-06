// Place all the behaviors and hooks related to the matching controller here.
// All this logic will automatically be available in application.js.

function rpubs_showLogin() {
  $('#login_message').hide();
  $('#login').modal();
}
function rpubs_logout() {
  $.ajax('/auth/logout.json', {
    type: 'POST',
    dataType: 'json'})
    .done(function(data) { window.location.reload(); })
    .fail(function(data) { debugger; });
}
function rpubs_login(username, password) {
  $('#login_message').hide();
  $.ajax('/auth/login.text', {
           type: 'POST',
           dataType: 'text',
           data: {
             username: username,
             password: password
           }})
    .done(function(data) { window.location.reload(); })
    .fail(function(jqXHR, textStatus, errorThrown) {
      $('#login_message').show();
      $('#login_message').text(jqXHR.responseText);
    });
}
$(function() {
  var doLogin = function() {
    var username = $('#login_username').val();
    var password = $('#login_password').val();
    $('#login_password').val('');
    rpubs_login(username, password);
  };

  $('#login').on('hidden', function() {
    $('#login_password').val('');
  });
  $('#login').on('shown', function() {
    $('#login_username').focus();
  });
  $('#login-modal-submit').click(doLogin);
  $('#login-modal-cancel').click(function() {
    $('#login_password').val('');
    $('#login').modal('hide');
  });

  var keyHandler = function(event) {
    if (event.which == 13)
      doLogin();
  };
  $('#login_username').keydown(keyHandler);
  $('#login_password').keydown(keyHandler);
});
