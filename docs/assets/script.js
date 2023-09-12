(function($) {
    $('.rst-content a[href^="http"]').attr('target', '_blank');

    $(`.rst-content #faq ~ details${window.location.hash}`).attr('open', true);
})(jQuery);
